from __future__ import annotations

import csv
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import basix.ufl
import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem


def parse_vtu_mesh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Читает узлы и тетраэдры из ASCII VTU-файла."""
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    if piece is None:
      raise RuntimeError(f"VTU mesh piece not found in {path}")

    points_text = piece.find("./Points/DataArray")
    if points_text is None or points_text.text is None:
      raise RuntimeError(f"VTU points not found in {path}")

    cells_text = piece.find("./Cells/DataArray[@Name='connectivity']")
    if cells_text is None or cells_text.text is None:
      raise RuntimeError(f"VTU connectivity not found in {path}")

    points = np.fromstring(points_text.text, sep=" ", dtype=np.float64).reshape(-1, 3)
    cells = np.fromstring(cells_text.text, sep=" ", dtype=np.int64).reshape(-1, 4)
    return points, cells


def read_frame_loads(path: Path) -> list[dict[str, float]]:
    """Читает по кадрам нагрузки, подготовленные C++."""
    with path.open("r", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        frames = []
        for row in reader:
            frame = {key: float(value) for key, value in row.items() if key not in {"frame_id"}}
            frame["frame_id"] = int(row["frame_id"])
            frames.append(frame)
    return frames


def write_box_vtu(
    path: Path,
    points: np.ndarray,
    cells: np.ndarray,
    density: float,
    wall_thickness: float,
    box_rgb: tuple[int, int, int],
    speed: float,
    roof_deflection: float,
    roof_center_x: float,
    roof_center_y: float,
    floor_deflection: float,
    side_deflection: float,
):
    """Записывает деформированный бокс в VTU-файл."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        out.write("  <UnstructuredGrid>\n")
        out.write(f'    <Piece NumberOfPoints="{len(points)}" NumberOfCells="{len(cells)}">\n')
        out.write("      <Points>\n")
        out.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for point in points:
            out.write(f"          {point[0]:.6f} {point[1]:.6f} {point[2]:.6f}\n")
        out.write("        </DataArray>\n")
        out.write("      </Points>\n")

        out.write('      <CellData Scalars="wall_thickness">\n')
        for name, value in [
            ("wall_thickness", wall_thickness),
            ("density", density),
            ("mass", 0.0),
            ("speed", speed),
            ("roof_deflection", roof_deflection),
            ("roof_center_x", roof_center_x),
            ("roof_center_y", roof_center_y),
            ("floor_deflection", floor_deflection),
            ("side_deflection", side_deflection),
        ]:
            out.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
            for _ in range(len(cells)):
                out.write(f"          {value:.6f}\n")
            out.write("        </DataArray>\n")

        out.write('        <DataArray type="UInt8" Name="box_rgb" NumberOfComponents="3" format="ascii">\n')
        for _ in range(len(cells)):
            out.write(f"          {box_rgb[0]} {box_rgb[1]} {box_rgb[2]}\n")
        out.write("        </DataArray>\n")
        out.write("      </CellData>\n")

        out.write("      <Cells>\n")
        out.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for tet in cells:
            out.write(f"          {tet[0]} {tet[1]} {tet[2]} {tet[3]}\n")
        out.write("        </DataArray>\n")

        out.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        for i in range(1, len(cells) + 1):
            out.write(f"          {4 * i}\n")
        out.write("        </DataArray>\n")

        out.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in range(len(cells)):
            out.write("          10\n")
        out.write("        </DataArray>\n")
        out.write("      </Cells>\n")
        out.write("    </Piece>\n")
        out.write("  </UnstructuredGrid>\n")
        out.write("</VTKFile>\n")


def write_box_pvd(path: Path, frames: list[dict[str, float]]):
    """Собирает все кадры бокса в одну серию для ParaView."""
    with path.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0"?>\n')
        out.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        out.write("  <Collection>\n")
        for frame in frames:
            out.write(
                f'    <DataSet timestep="{frame["time"]:.6f}" group="" part="0" '
                f'file="box/box_{frame["frame_id"]:04d}.vtu"/>\n'
            )
        out.write("  </Collection>\n")
        out.write("</VTKFile>\n")


def compute_deformation_metrics(
    reference_points: np.ndarray,
    deformed_points: np.ndarray,
    min_point: np.ndarray,
    max_point: np.ndarray,
) -> dict[str, float]:
    """Оценивает смятие каждой грани по узлам, лежащим близко к этой грани."""
    box_size = max_point - min_point
    tolerance = max(0.01, 0.015 * float(np.min(box_size)))

    roof_mask = np.isclose(reference_points[:, 2], max_point[2], atol=tolerance)
    floor_mask = np.isclose(reference_points[:, 2], min_point[2], atol=tolerance)
    xpos_mask = np.isclose(reference_points[:, 0], max_point[0], atol=tolerance)
    xneg_mask = np.isclose(reference_points[:, 0], min_point[0], atol=tolerance)
    ypos_mask = np.isclose(reference_points[:, 1], max_point[1], atol=tolerance)
    yneg_mask = np.isclose(reference_points[:, 1], min_point[1], atol=tolerance)

    metrics = {
        "roof": 0.0,
        "floor": 0.0,
        "x_pos": 0.0,
        "x_neg": 0.0,
        "y_pos": 0.0,
        "y_neg": 0.0,
    }

    if np.any(roof_mask):
        metrics["roof"] = max(
            0.0,
            float(np.max(reference_points[roof_mask, 2] - deformed_points[roof_mask, 2])),
        )

    if np.any(floor_mask):
        metrics["floor"] = max(
            0.0,
            float(np.max(deformed_points[floor_mask, 2] - reference_points[floor_mask, 2])),
        )

    if np.any(xpos_mask):
        metrics["x_pos"] = max(
            0.0,
            float(np.max(reference_points[xpos_mask, 0] - deformed_points[xpos_mask, 0])),
        )
    if np.any(xneg_mask):
        metrics["x_neg"] = max(
            0.0,
            float(np.max(deformed_points[xneg_mask, 0] - reference_points[xneg_mask, 0])),
        )
    if np.any(ypos_mask):
        metrics["y_pos"] = max(
            0.0,
            float(np.max(reference_points[ypos_mask, 1] - deformed_points[ypos_mask, 1])),
        )
    if np.any(yneg_mask):
        metrics["y_neg"] = max(
            0.0,
            float(np.max(deformed_points[yneg_mask, 1] - reference_points[yneg_mask, 1])),
        )

    return metrics


def update_crush_state(
    crush_state: dict[str, float],
    metrics: dict[str, float],
    frame: dict[str, float],
):
    """Обновляет накопленное смятие, которое уже не должно исчезать в следующих кадрах."""
    roof_load = float(frame["roof_load"])
    side_load = max(float(frame["side_x_load"]), float(frame["side_y_load"]))

    # Крыша запоминает самый сильный прогиб и дальше уже не распрямляется.
    crush_state["roof_sink"] = max(
        crush_state["roof_sink"],
        0.92 * metrics["roof"],
    )

    # Боковые стенки завязываем не только на прямой боковой удар,
    # но и на смятие крыши, чтобы все 4 грани втягивались более согласованно.
    target_side = max(
        metrics["x_pos"],
        metrics["x_neg"],
        metrics["y_pos"],
        metrics["y_neg"],
        0.22 * crush_state["roof_sink"],
        0.00002 * roof_load,
        0.00012 * side_load,
    )
    crush_state["side_sink"] = max(crush_state["side_sink"], 0.88 * target_side)

    # Днище тоже может слегка изменить форму из-за опоры на землю.
    crush_state["floor_raise"] = max(crush_state["floor_raise"], 0.55 * metrics["floor"])

    # Запоминаем последнее значимое место удара по крыше.
    if roof_load > 1e-6:
        crush_state["roof_center_x"] = 0.80 * crush_state["roof_center_x"] + 0.20 * float(frame["roof_center_x"])
        crush_state["roof_center_y"] = 0.80 * crush_state["roof_center_y"] + 0.20 * float(frame["roof_center_y"])


def build_permanent_crush_field(
    reference_points: np.ndarray,
    min_point: np.ndarray,
    max_point: np.ndarray,
    crush_state: dict[str, float],
) -> np.ndarray:
    """Строит постоянную геометрию смятого кузова из накопленных параметров."""
    center = 0.5 * (min_point + max_point)
    half = 0.5 * (max_point - min_point)
    local = reference_points - center

    displacement = np.zeros_like(reference_points)

    # Нормированные координаты по сторонам нужны как плавные веса деформации.
    x_ratio = np.clip(np.abs(local[:, 0]) / max(half[0], 1e-9), 0.0, 1.0)
    y_ratio = np.clip(np.abs(local[:, 1]) / max(half[1], 1e-9), 0.0, 1.0)
    z_upper = np.clip(local[:, 2] / max(half[2], 1e-9), 0.0, 1.0)
    z_lower = np.clip(-local[:, 2] / max(half[2], 1e-9), 0.0, 1.0)

    roof_patch_x = max(0.45, 0.28 * float(max_point[0] - min_point[0]))
    roof_patch_y = max(0.30, 0.32 * float(max_point[1] - min_point[1]))
    dx = local[:, 0] - crush_state["roof_center_x"]
    dy = local[:, 1] - crush_state["roof_center_y"]
    roof_patch = np.exp(
        -(dx * dx) / (2.0 * roof_patch_x * roof_patch_x)
        -(dy * dy) / (2.0 * roof_patch_y * roof_patch_y)
    )

    # Крыша проминается и локально в точке удара, и более широко по всей верхней части.
    roof_sink = crush_state["roof_sink"]
    displacement[:, 2] -= roof_sink * (
        0.55 * z_upper * z_upper + 0.45 * z_upper * z_upper * roof_patch
    )

    # Все 4 боковые стенки дополнительно сминаются внутрь.
    side_sink = crush_state["side_sink"]
    side_weight_x = x_ratio * x_ratio * (0.30 + 0.70 * z_upper)
    side_weight_y = y_ratio * y_ratio * (0.30 + 0.70 * z_upper)
    displacement[:, 0] -= np.sign(local[:, 0]) * side_sink * side_weight_x
    displacement[:, 1] -= np.sign(local[:, 1]) * side_sink * side_weight_y

    # Когда крыша уходит вниз, верхняя часть стенок тоже оседает.
    wall_drop = 0.42 * roof_sink
    wall_weight = np.maximum(side_weight_x, side_weight_y)
    displacement[:, 2] -= wall_drop * wall_weight

    # Днище меняется слабее, потому что основное смятие идёт сверху.
    floor_raise = crush_state["floor_raise"]
    displacement[:, 2] += floor_raise * z_lower * z_lower

    # Ограничиваем диапазон остаточной формы, чтобы геометрия не становилась совсем неустойчивой.
    displacement[:, 0] = np.clip(displacement[:, 0], -0.95, 0.95)
    displacement[:, 1] = np.clip(displacement[:, 1], -0.95, 0.95)
    displacement[:, 2] = np.clip(displacement[:, 2], -1.70, 0.15)

    return displacement


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: fenics_box_solver.py <results_dir>", file=sys.stderr)
        return 1

    results_dir = Path(sys.argv[1]).resolve()
    mesh_path = results_dir / "mesh" / "box_shell.vtu"
    loads_path = results_dir / "box_loads.csv"
    output_box_dir = results_dir / "box"

    points, cells = parse_vtu_mesh(mesh_path)
    frames = read_frame_loads(loads_path)
    output_box_dir.mkdir(parents=True, exist_ok=True)

    # Создаём сетку dolfinx прямо из уже готовых узлов и тетраэдров.
    coord_element = basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))
    domain = mesh.create_mesh(MPI.COMM_SELF, cells, coord_element, points)

    # Пространство векторных перемещений.
    V = fem.functionspace(domain, ("Lagrange", 1, (3,)))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Полноценной пластичности и потери устойчивости пока нет.
    # Поэтому используем эффективную жёсткость ниже модуля стали,
    # чтобы грубо имитировать смятие тонкостенного кузова после потери устойчивости.
    young_modulus = 1.5e9
    poisson_ratio = 0.30
    mu = young_modulus / (2.0 * (1.0 + poisson_ratio))
    lmbda = (
        young_modulus
        * poisson_ratio
        / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
    )

    def epsilon(w):
        return ufl.sym(ufl.grad(w))

    def sigma(w):
        return 2.0 * mu * epsilon(w) + lmbda * ufl.tr(epsilon(w)) * ufl.Identity(3)

    # Геометрические размеры нужны для пятна нагрузки и граничных условий.
    min_point = points.min(axis=0)
    max_point = points.max(axis=0)
    box_size = max_point - min_point
    z_min = float(min_point[2])
    z_max = float(max_point[2])

    roof_sigma_x = max(0.35, 0.20 * float(box_size[0]))
    roof_sigma_y = max(0.25, 0.24 * float(box_size[1]))
    side_sigma_z = max(0.35, 0.35 * float(box_size[2]))
    top_area = max(float(box_size[0] * box_size[1]), 1e-8)
    side_area_x = max(float(box_size[1] * box_size[2]), 1e-8)
    side_area_y = max(float(box_size[0] * box_size[2]), 1e-8)

    # Эти коэффициенты переводят короткий удар воды в более заметное смятие упрощённого кузова.
    impact_force_scale = 30.0
    roof_uniform_fraction = 0.55
    roof_patch_fraction = 0.45
    wall_crush_fraction = 0.16

    # Постоянные будут меняться от кадра к кадру без перекомпиляции формы.
    roof_load = fem.Constant(domain, 0.0)
    roof_center_x = fem.Constant(domain, 0.0)
    roof_center_y = fem.Constant(domain, 0.0)
    side_x_load = fem.Constant(domain, 0.0)
    side_y_load = fem.Constant(domain, 0.0)

    x = ufl.SpatialCoordinate(domain)

    # Мягкое распределение удара по крыше.
    roof_area = max(2.0 * math.pi * roof_sigma_x * roof_sigma_y, 1e-8)
    roof_uniform_pressure = roof_uniform_fraction * roof_load / top_area
    roof_patch_pressure = roof_patch_fraction * roof_load / roof_area
    roof_patch = ufl.exp(
        -((x[0] - roof_center_x) ** 2) / (2.0 * roof_sigma_x * roof_sigma_x)
        -((x[1] - roof_center_y) ** 2) / (2.0 * roof_sigma_y * roof_sigma_y)
    )
    traction_top = ufl.as_vector(
        (0.0, 0.0, -(roof_uniform_pressure + roof_patch_pressure * roof_patch))
    )

    # Боковые стенки тоже начинают сминаться, когда крыша уходит вниз.
    side_pressure_x = side_x_load / side_area_x + wall_crush_fraction * roof_load / side_area_x
    side_pressure_y = side_y_load / side_area_y + wall_crush_fraction * roof_load / side_area_y
    side_profile_z = 0.55 + 0.45 * ufl.exp(
        -((x[2] - z_max) ** 2) / (2.0 * side_sigma_z * side_sigma_z)
    )
    traction_x_pos = ufl.as_vector((-side_pressure_x * side_profile_z, 0.0, 0.0))
    traction_x_neg = ufl.as_vector((side_pressure_x * side_profile_z, 0.0, 0.0))
    traction_y_pos = ufl.as_vector((0.0, -side_pressure_y * side_profile_z, 0.0))
    traction_y_neg = ufl.as_vector((0.0, side_pressure_y * side_profile_z, 0.0))

    # Нижнюю грань считаем зафиксированной землёй
    fdim = domain.topology.dim - 1
    bottom_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda X: np.isclose(X[2], z_min)
    )
    bottom_dofs = fem.locate_dofs_topological(V, fdim, bottom_facets)
    bc = fem.dirichletbc(np.array((0.0, 0.0, 0.0), dtype=np.float64), bottom_dofs, V)

    top_facets = mesh.locate_entities_boundary(domain, fdim, lambda X: np.isclose(X[2], z_max))
    xpos_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda X: np.isclose(X[0], max_point[0])
    )
    xneg_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda X: np.isclose(X[0], min_point[0])
    )
    ypos_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda X: np.isclose(X[1], max_point[1])
    )
    yneg_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda X: np.isclose(X[1], min_point[1])
    )

    tagged_facets = np.hstack(
        [top_facets, xpos_facets, xneg_facets, ypos_facets, yneg_facets]
    ).astype(np.int32)
    tagged_values = np.hstack(
        [
            np.full(len(top_facets), 1, dtype=np.int32),
            np.full(len(xpos_facets), 2, dtype=np.int32),
            np.full(len(xneg_facets), 3, dtype=np.int32),
            np.full(len(ypos_facets), 4, dtype=np.int32),
            np.full(len(yneg_facets), 5, dtype=np.int32),
        ]
    )
    order = np.argsort(tagged_facets)
    facet_tags = mesh.meshtags(domain, fdim, tagged_facets[order], tagged_values[order])
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

    # Базовая форма линейной упругости
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
    L = (
        ufl.dot(traction_top, v) * ds(1)
        + ufl.dot(traction_x_pos, v) * ds(2)
        + ufl.dot(traction_x_neg, v) * ds(3)
        + ufl.dot(traction_y_pos, v) * ds(4)
        + ufl.dot(traction_y_neg, v) * ds(5)
    )

    problem = LinearProblem(
        a,
        L,
        bcs=[bc],
        petsc_options_prefix="box_elasticity_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )

    # Цвет бокса меняем на стальной серый, чтобы не сливался с землёй.
    box_rgb = (160, 170, 182)
    wall_thickness = 0.0008

    # Здесь храним накопленное смятие кузова
    crush_state = {
        "roof_sink": 0.0,
        "side_sink": 0.0,
        "floor_raise": 0.0,
        "roof_center_x": 0.0,
        "roof_center_y": 0.0,
    }

    # Упругая часть тоже остаётся видимой, но уже слабее остаточной
    elastic_visible_fraction = 0.10

    for frame in frames:
        # Эти значения приходят из C++ как эквивалентные нагрузки от воды.
        roof_load.value = impact_force_scale * float(frame["roof_load"])
        roof_center_x.value = float(frame["roof_center_x"])
        roof_center_y.value = float(frame["roof_center_y"])
        side_x_load.value = impact_force_scale * float(frame["side_x_load"])
        side_y_load.value = impact_force_scale * float(frame["side_y_load"])

        # Решаем упругую задачу для текущего кадра.
        solution = problem.solve()
        elastic_displacement = solution.x.array.reshape(-1, 3)

        # Сдвиг всего бокса оставляем минимальным
        shift = np.array(
            [
                0.0,
                0.0,
                min(0.0, float(frame["tz"])),
            ],
            dtype=np.float64,
        )

        # Сначала оцениваем чистый упругий отклик, а затем обновляем
        # накопленное смятие, которое уже не должно исчезать дальше.
        elastic_points = points + elastic_displacement + shift
        elastic_metrics = compute_deformation_metrics(points, elastic_points, min_point, max_point)
        update_crush_state(crush_state, elastic_metrics, frame)

        permanent_displacement = build_permanent_crush_field(points, min_point, max_point, crush_state)
        total_displacement = permanent_displacement + elastic_visible_fraction * elastic_displacement
        deformed_points = points + total_displacement + shift

        # Земля допускает только небольшое вдавливание корпуса вниз.
        # Это убирает неестественный уход сетки глубоко под плоскость опоры.
        max_ground_sink = 0.03
        deformed_points[:, 2] = np.maximum(deformed_points[:, 2], z_min - max_ground_sink)

        final_metrics = compute_deformation_metrics(points, deformed_points, min_point, max_point)
        roof_deflection = final_metrics["roof"]
        floor_deflection = final_metrics["floor"]
        side_deflection = max(
            final_metrics["x_pos"],
            final_metrics["x_neg"],
            final_metrics["y_pos"],
            final_metrics["y_neg"],
        )

        write_box_vtu(
            output_box_dir / f"box_{frame['frame_id']:04d}.vtu",
            deformed_points,
            cells,
            density=7800.0,
            wall_thickness=wall_thickness,
            box_rgb=box_rgb,
            speed=0.0,
            roof_deflection=roof_deflection,
            roof_center_x=crush_state["roof_center_x"],
            roof_center_y=crush_state["roof_center_y"],
            floor_deflection=floor_deflection,
            side_deflection=side_deflection,
        )

    write_box_pvd(results_dir / "box_series.pvd", frames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
