/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//

#include "./TriangleList.h"
#include "./tricase.cxx"

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>





// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
  //  return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
     idx[0] = pointId%dims[0];
     idx[1] = (pointId/dims[0])%dims[1];
     idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
     idx[0] = cellId%(dims[0]-1);
     idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
     idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}


class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    cerr<<"Testing tricase"<<triCase[0][1]<<endl;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    //list of edges per cell

    static int edges[12][2] = {{0,1}, {1,3}, {2,3}, {0,2}, {4,5}, {5,7},
				{6,7}, {4,6}, {0,4}, {1,5}, {2,6}, {3,7}};

    //initalizing values for the vertex,logical indices x y and z, and the 1 or 0 to identify the case

    int b0, b1, b2, b3, b4, b5, b6, b7;
    int b[8];
    int b_logical[8][3];

   //The isovalue for proj6b is 3.2

    float iso = 3.2;
    int num_cells = GetNumberOfCells(dims);
    
    //now I need to add all triangles with the isovalue 3.2 to my list tl
    //first taking each cell from 0->dims[0]-1*dims[1]-1

    //Everything up to this point seens to work!

    TriangleList tl;
    int n;
    for (int i=0;i<1000;i++){
	cerr<<"CELL "<<i<<endl;	

	int cell_idx[3];

	//Get logical cell index for each cell

	GetLogicalCellIndex(cell_idx ,i, dims);


	//conventions per the lecture slides shows that 
	//b0->b1 is along the x axis
	//b0->b2 is along the y axis
	//b0->b4 is along the z axis
	//now get each corner starting w lower left

	b_logical[0][0] = cell_idx[0];
	b_logical[0][1] = cell_idx[1];
	b_logical[0][2] =  cell_idx[2];
	b[0] = GetPointIndex(b_logical[0], dims);
	
	b_logical[1][0] = cell_idx[0]+1;
	b_logical[1][1] = cell_idx[1];
	b_logical[1][2] = cell_idx[2];
	b[1] = GetPointIndex(b_logical[1], dims);

	b_logical[2][0] = cell_idx[0];
	b_logical[2][1] = cell_idx[1]+1;
	b_logical[2][2] = cell_idx[2];
	b[2] = GetPointIndex(b_logical[2], dims);
	
	b_logical[3][0] = cell_idx[0]+1;
	b_logical[3][1] = cell_idx[1]+1;
	b_logical[3][2] = cell_idx[2];
	b[3] = GetPointIndex(b_logical[3], dims);

	b_logical[4][0] = cell_idx[0];
	b_logical[4][1] = cell_idx[1]; 
	b_logical[4][2] = cell_idx[2]+1;
	b[4] = GetPointIndex(b_logical[4], dims);

	b_logical[5][0] = cell_idx[0]+1;
	b_logical[5][1] = cell_idx[1];
	b_logical[5][2] =  cell_idx[2]+1;
	b[5] = GetPointIndex(b_logical[5], dims);

	b_logical[6][0] = cell_idx[0];
	b_logical[6][1] =  cell_idx[1]+1;
	b_logical[6][2] = cell_idx[2]+1;
	b[6] = GetPointIndex(b_logical[6], dims);

	b_logical[7][0] = cell_idx[0]+1; 
	b_logical[7][1] = cell_idx[1]+1;
	b_logical[7][2] = cell_idx[2]+1;
	b[7] = GetPointIndex(b_logical[7], dims);

	//get value for F(v) where v = each vertice

	float A = F[b[0]];
	float B = F[b[1]];
	float C = F[b[2]];
	float D = F[b[3]];
	float E = F[b[4]];
	float G = F[b[5]];
	float H = F[b[6]];
	float I = F[b[7]];
//	cerr<<A<<" "<<B << " "<<C<<" "<<D<<"ABCD-EGHI "<<E<<" "<<G<<" "<<H<<" "<<I<<endl;

	//identify the cases 0-255 using binary.
	//F(b0) = A

	if (A >= iso)
		b0 = 1;
	else{ b0 = 0; }
	
	if (B >= iso)
		b1 = 1;
	else{ b1 = 0; }
	
	if (C >= iso)
		b2 = 1;
	else{ b2 = 0; }
	
	if (D >= iso)
		b3 = 1;
	else{ b3 = 0; }
	
	if (E >= iso)
		b4 = 1;
	else{ b4 = 0; }
	
	if (G >= iso)
		b5 = 1;
	else{ b5 = 0; }
	
	if (H >= iso)
		b6 = 1;
	else{ b6 = 0; }
	
	if (I >= iso)
		b7 = 1;
	else{ b7 = 0; }

	//using the binary number to convert to a case number
	
	int icase = (1*b0)+(2*b1)+(4*b2)+(8*b3)+(16*b4)+(32*b5)+(64*b6)+(128*b7);	
        j=0;
//	cerr<<b0<<" "<<b1 << " "<<b2<<" "<<b3<<"0123-4567 "<<b4<<" "<<b5<<" "<<b6<<" "<<b7<<" "<<icase <<endl;
	
	//now for each cell case while triCase[icase][i] != -1 
	while (triCase[icase][j] != -1)
	
	    //AddTriangle(X1,Y1,Z1, X2, Y2, Z2, X3, Y3, Z3);
	    {
	
	    float X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3;
	
	    //parse 3 edges at a time, if we encounter -1 we are done

	    int edge1[2],edge2[2],edge3[2];
	    int e1, e2, e3;
	    float e1_t, e2_t, e3_t;
	    e1 = triCase[icase][j];
	    e2 = triCase[icase][j+1];
	    e3 = triCase[icase][j+2];

	    //i have an edge e1-3 which correspond to edges of a cube
	    //edges[][] corresponds to the connecting vertices of this edge

	    edge1[0] = edges[e1][0];
	    edge1[1] = edges[e1][1];
	    edge2[0] = edges[e2][0];
	    edge2[1] = edges[e2][1];
	    edge3[0] = edges[e3][0];
	    edge3[1] = edges[e3][1];
	    //these need to be equal to A or B
	    float E1F1 = F[b[edge1[0]]];
	    float E1F2 = F[b[edge1[1]]];
	    float E2F1 = F[b[edge2[0]]];
	    float E2F2 = F[b[edge2[1]]];
	    float E3F1 = F[b[edge3[0]]];
	    float E3F2 = F[b[edge3[1]]];
		
	   
	    //now i want to interpolate the points between each edge
	    //say we have case 1, so 0, 3, 8
	    //that gives triangle of vertices (0,1),(0,2),(0,4)
	    //so we need the logical indices of the points 0,1
	    
	 

	     e1_t =  (iso - E1F1) / (E1F2 - E1F1);
	    cerr<<"("<< iso<<"-"<<E1F1<<" )"<<"/ "<< "(" <<E1F2<<"-"<<E1F1<<endl;
	     e2_t =  (iso - E2F1) / (E2F2 - E2F1);
	    cerr<<"("<< iso<<"-"<<E2F1<<" )"<<"/ "<< "(" <<E2F2<<"-"<<E2F1<<endl;

	     e3_t =  (iso - E3F1) / (E3F2 - E3F1);
	    cerr<<"("<< iso<<"-"<<E3F1<<" )"<<"/ "<< "(" <<E3F2<<"-"<<E3F1<<endl;

//            cerr<<e1_t<<" "<<e2_t<< " "<<e3_t<< endl;

//	    cerr<< e1<<" "<< e1_t<<" F->" <<F[edge1[0]] << " "<< F[edge1[1]]<<endl;
	    //edge on z axis? if so interpolate on the z axis and get scale the other coordinates

	    if (e1 >= 8 && e1 < 12)
		{
		//I want x1 to be the first x coordinate on the edges[0] vertex if the edge is 1
	        e1_t = Z[b_logical[edges[e1][0]][2]] + e1_t*(Z[b_logical[edges[e1][1]][2]] - Z[b_logical[edges[e1][0]][2]]);

		X1 = X[b_logical[edges[e1][0]][0]];
		Y1 = Y[b_logical[edges[e1][0]][1]];
		Z1 = e1_t;
		}
	    //x axis?
	    else if (e1%2 == 0)
		{
		e1_t = X[b_logical[edges[e1][0]][0]] + e1_t*(X[b_logical[edges[e1][1]][0]] - X[b_logical[edges[e1][0]][0]]);

		X1 = e1_t;
		Y1 = Y[b_logical[edges[e1][0]][1]];
		Z1 = Z[b_logical[edges[e1][0]][2]];
		}
	    //y axis?
	    else if (e1%2 == 1)
		{
	        e1_t = Y[b_logical[edges[e1][0]][2]] + e1_t*(Y[b_logical[edges[e1][1]][2]] - Y[b_logical[edges[e1][0]][2]]);

		X1 = X[b_logical[edges[e1][0]][0]];
		Y1 = e1_t;
		Z1 = Z[b_logical[edges[e1][0]][2]];
		}
	   	


	    if (e2 >= 8 && e2 < 12)
		{
		e2_t = Z[b_logical[edges[e2][0]][2]] + e2_t*(Z[b_logical[edges[e2][1]][2]] - Z[b_logical[edges[e2][0]][2]]);

		X2 = X[b_logical[edges[e2][0]][0]];
		Y2 = Y[b_logical[edges[e2][0]][1]];
		Z2 = e2_t;
		}
	    else if (e2%2 == 0)
		{
	        e2_t = X[b_logical[edges[e2][0]][2]] + e2_t*(X[b_logical[edges[e2][1]][2]] - X[b_logical[edges[e2][0]][2]]);

		X2 = e2_t;
		Y2 = Y[b_logical[edges[e2][0]][1]];
		Z2 = Z[b_logical[edges[e2][0]][2]];
		}
	    else if (e2%2 == 1)
		{
	        e2_t = Y[b_logical[edges[e2][0]][2]] + e2_t*(Y[b_logical[edges[e2][1]][2]] - Y[b_logical[edges[e2][0]][2]]);

		X2 = X[b_logical[edges[e2][0]][0]];
		Y2 = e2_t;
		Z2 = Z[b_logical[edges[e2][0]][2]];
		}

	    
	    if (e3 >= 8 && e3 <12)
		{
	        e3_t = Z[b_logical[edges[e3][0]][2]] + e3_t*(Z[b_logical[edges[e3][1]][2]] - Z[b_logical[edges[e3][0]][2]]);

		X3 = X[b_logical[edges[e3][0]][0]];
		Y3 = Y[b_logical[edges[e3][0]][1]];
		Z3 = e3_t;		
		}
	    else if (e3%2 == 0)
		{
	        e3_t = X[b_logical[edges[e3][0]][2]] + e3_t*(X[b_logical[edges[e3][1]][2]] - X[b_logical[edges[e3][0]][2]]);

		X3 = e3_t;
		Y3 = Y[b_logical[edges[e3][0]][1]];
		Z3 = Z[b_logical[edges[e3][0]][2]];
		
		}
	   else if (e3%2 == 1)
		{
	        e3_t = Y[b_logical[edges[e3][0]][2]] + e3_t*(Y[b_logical[edges[e3][1]][2]] - Y[b_logical[edges[e3][0]][2]]);

		X3 = X[b_logical[edges[e3][0]][0]];
		Y3 = e3_t;
		Z3 = Z[b_logical[edges[e3][0]][2]];
		
		}
	    tl.AddTriangle(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3);

	    j+=3;
	
            }//endwhile


	}//endfor1
    cerr<<"Finshed for # of cells"<< (dims[0]-1)*(dims[1]-1)*(dims[2]-1)<<endl;
//    cerr<< triCase[0][0]<<endl;
    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */
    vtkPolyData *pd = vtkPolyData::New();
    pd = tl.MakePolyData();
    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}
