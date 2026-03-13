#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

class CalcNode
{
friend class CalcMesh;

protected:
    double x;
    double y;
    double z;

    double smth;

    double vx;
    double vy;
    double vz;

public:
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0)
    {
    }

    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
            : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz)
    {
    }

    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};

class Element
{
friend class CalcMesh;

protected:
    unsigned long nodesIds[4];
};

class CalcMesh
{
protected:
    vector<CalcNode> points;
    vector<Element> elements;

public:
    CalcMesh(const vector<double>& nodesCoords, const vector<size_t>& tetrsPoints) {
        points.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];

            double smth = sqrt(pointX*pointX + pointY*pointY + pointZ*pointZ);

            double stretchX = 0.05;
            double stretchYZ = -stretchX / 2;

            double vx = stretchX * pointX;
            double vy = stretchYZ * pointY;
            double vz = stretchYZ * pointZ;

            points[i] = CalcNode(pointX, pointY, pointZ, smth, vx, vy, vz);
        }

        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    void doTimeStep(double tau) {
        for(unsigned int i = 0; i < points.size(); i++) {
            points[i].move(tau);
        }
    }

    void snapshot(unsigned int snap_number) {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        unsigned int number = (unsigned int)points.size();
        for(unsigned int i = 0; i < number; i++) {
            dumpPoints->InsertNextPoint(points[i].x, points[i].y, points[i].z);

            double _vel[3] = {points[i].vx, points[i].vy, points[i].vz};
            vel->InsertNextTuple(_vel);

            smth->InsertNextValue(points[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);

        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        string fileName = "dynamic_stl_obj_3D" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("dynamic_stl_obj_3D");

    try {
        gmsh::merge("t9.msh"); 
    } catch(...) {
        gmsh::logger::write("error!!!!!!!!");
        gmsh::finalize();
        return -1;
    }

    vector<double> nodesCoord;
    vector<size_t> nodeTags;
    vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    vector<int> elementTypes;
    vector<vector<size_t>> elementTags;
    vector<vector<size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    vector<size_t>* tetrsNodesTags = nullptr;
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] == GMSH_TETR_CODE) {
            tetrsNodesTags = &elementNodeTags[i];
            break;
        }
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    double tau = 0.05;

    mesh.snapshot(0);

    for(unsigned int step = 1; step <= 100; step++) {
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    return 0;
}