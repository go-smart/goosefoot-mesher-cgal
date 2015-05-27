/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOFFWriter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOFFWriter.h"

#include "vtkByteSwap.h"
#include "vtkCellArray.h"
#include "vtkErrorCode.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include <vtksys/SystemTools.hxx>

#if !defined(_WIN32) || defined(__CYGWIN__)
# include <unistd.h> /* unlink */
#else
# include <io.h> /* unlink */
#endif

vtkStandardNewMacro(vtkOFFWriter);

static char header[]="OFF";

vtkOFFWriter::vtkOFFWriter()
{
  this->FileType = VTK_ASCII;
  strcpy(this->Header, header);
}

void vtkOFFWriter::WriteData()
{
  vtkPoints *pts;
  vtkCellArray *polys;
  vtkPolyData *input = this->GetInput();

  polys = input->GetPolys();
  pts = input->GetPoints();
  if (pts == NULL || polys == NULL )
    {
    vtkErrorMacro(<<"No data to write!");
    return;
    }

  if ( this->FileName == NULL)
    {
    vtkErrorMacro(<< "Please specify FileName to write");
    this->SetErrorCode(vtkErrorCode::NoFileNameError);
    return;
    }

  if ( this->FileType == VTK_BINARY )
    {
    this->WriteBinaryOFF(pts,polys);
    if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
      {
      vtkErrorMacro("Ran out of disk space; deleting file: "
                    << this->FileName);
      unlink(this->FileName);
      }
    }
  else
    {
    this->WriteAsciiOFF(pts,polys);
    if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError)
      {
      vtkErrorMacro("Ran out of disk space; deleting file: "
                    << this->FileName);
      unlink(this->FileName);
      }
    }
}

void vtkOFFWriter::WriteAsciiOFF(vtkPoints *pts, vtkCellArray *polys)
{
  FILE *fp;
  double n[3], v1[3], v2[3], v3[3];
  vtkIdType npts = 0, ncells = 0;
  vtkIdType *indx = 0;
  
  if ((fp = fopen(this->FileName, "w")) == NULL)
    {
    vtkErrorMacro(<< "Couldn't open file: " << this->FileName);
    this->SetErrorCode(vtkErrorCode::CannotOpenFileError);
    return;
    }
//
//  Write header
//
  vtkDebugMacro("Writing ASCII off file");

  npts   = pts->GetNumberOfPoints();
  ncells = polys->GetNumberOfCells();

  if (fprintf (fp, "OFF\n%d %d 0\n", npts, ncells) < 0)
    {
    fclose(fp);
    this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
    return;
    }

  for (*indx = 0; pts->GetPoint(*indx, v1); *indx < npts )
    {
    if (fprintf (fp, "%.6g %.6g %.6g\n", v1[0], v1[1], v2[2]) < 0)
      {
      fclose(fp);
      this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
      return;
      }

  for (polys->InitTraversal(); polys->GetNextCell(npts,indx); )
    {
    if (fprintf (fp, "3 %d %d %d\n",
                 indx[0], indx[1], indx[2]) < 0)
      {
      fclose(fp);
      this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
      return;
      }
    }

  fclose (fp);
}

//----------------------------------------------------------------------------
void vtkOFFWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
