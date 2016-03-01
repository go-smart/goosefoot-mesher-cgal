/**
 * mesher_cgal
 *
 * Copyright (C) 2013-  NUMA Engineering Ltd. (see AUTHORs file)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// http://www.vtk.org/Wiki/VTK/Examples/Cxx/StructuredGrid/StructuredGrid

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include "mesher_cgal.h"

using namespace mesherCGAL;

vtkSmartPointer<vtkStructuredGrid> structuredGrid;
vtkSmartPointer<vtkDoubleArray> allocationOrder;
vtkSmartPointer<vtkDoubleArray> characteristicLength;
vtkSmartPointer<vtkDoubleArray> needleDist;
vtkSmartPointer<vtkDoubleArray> location;

int mesherCGAL::visualization_set_allocation_order(int ix, int order, double cl, double needle_dist, double x, double y, double z) {
    allocationOrder->SetValue(ix, (double)order);
    characteristicLength->SetValue(ix, cl);
    needleDist->SetValue(ix, needle_dist);
    double xv[3] = {x, y, z};
    location->SetTuple(ix, xv);
	return EXIT_SUCCESS;
}

int mesherCGAL::visualization_save(std::string& filename) {
  // Write file
  vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(structuredGrid);
#else
  writer->SetInputData(structuredGrid);
#endif
  writer->Write();
  return EXIT_SUCCESS;
}

int mesherCGAL::visualization_create_structured_grid(double x0, double y0, double z0, int nx, int ny, int nz, double dx) {
  // Create a grid
  structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();
  allocationOrder =
      vtkSmartPointer<vtkDoubleArray>::New();
  needleDist =
      vtkSmartPointer<vtkDoubleArray>::New();
  characteristicLength =
      vtkSmartPointer<vtkDoubleArray>::New();
  location =
      vtkSmartPointer<vtkDoubleArray>::New();

  location->SetNumberOfComponents(3);
  location->SetName("Location");
  allocationOrder->SetNumberOfComponents(1);
  allocationOrder->SetName("Allocation Order");
  needleDist->SetNumberOfComponents(1);
  needleDist->SetName("Distance from Needle");
  characteristicLength->SetNumberOfComponents(1);
  characteristicLength->SetName("Characteristic Length");
  //allocationOrder->SetNumberOfTuples(nx * ny * nz);
  std::cout << nx << " " << ny << " " << nz << std::endl;
 
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

  double x, y, z;
 
  double zero[3] = {0., 0., 0.,};
  z = z0;
  for(unsigned int i = 0; i <= nz; i++)
    {
    z += dx;
    y = y0;
    for(unsigned int j = 0; j <= ny; j++)
      {
      y += dx;
      x = x0;
      for(unsigned int k = 0; k <= nx; k++)
        {
        x += dx;
        points->InsertNextPoint(x, y, z);
        if (i > 0 && j > 0 && k > 0) {
            allocationOrder->InsertNextValue(-1.0);
            location->InsertNextTuple(zero);
            characteristicLength->InsertNextValue(-1.0);
            needleDist->InsertNextValue(-1.0);
        }
        }
      }
    }
 
  // Specify the dimensions of the grid
  structuredGrid->SetDimensions(nx+1,ny+1,nz+1);
  structuredGrid->SetPoints(points);
  structuredGrid->GetCellData()->AddArray(allocationOrder);
  structuredGrid->GetCellData()->AddArray(location);
  structuredGrid->GetCellData()->AddArray(characteristicLength);
  structuredGrid->GetCellData()->AddArray(needleDist);

  // Create a mapper and actor
  vtkSmartPointer<vtkDataSetMapper> mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInputConnection(structuredGrid->GetProducerPort());
#else
  mapper->SetInputData(structuredGrid);
#endif
 
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
 
  return EXIT_SUCCESS;
}
