/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOFFWriter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOFFWriter - write stereo lithography files
// .SECTION Description
// vtkOFFWriter writes stereo lithography (.off) files in either ASCII or
// binary form. Stereo lithography files only contain triangles. If polygons
// with more than 3 vertices are present, only the first 3 vertices are
// written.  Use vtkTriangleFilter to convert polygons to triangles.

// .SECTION Caveats
// Binary files written on one system may not be readable on other systems.
// vtkOFFWriter uses VAX or PC byte ordering and swaps bytes on other systems.

#ifndef __vtkOFFWriter_h
#define __vtkOFFWriter_h

#include "vtkPolyDataWriter.h"

class VTK_IO_EXPORT vtkOFFWriter : public vtkPolyDataWriter
{
public:
  static vtkOFFWriter *New();
  vtkTypeMacro(vtkOFFWriter,vtkPolyDataWriter);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkOFFWriter();
  ~vtkOFFWriter() {};

  void WriteData();

  void WriteAsciiOFF(vtkPoints *pts, vtkCellArray *polys);
private:
  vtkOFFWriter(const vtkOFFWriter&);  // Not implemented.
  void operator=(const vtkOFFWriter&);  // Not implemented.
};

#endif

