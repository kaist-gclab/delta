/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the GPL license (see License.txt)
* ------------------------------------------------------------------------ */  

#include <stack>
#include <assert.h>
#include <vtkObjectFactory.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkIdListCollection.h>

#include "vtkSurfaceBase.h"

int vtkSurfaceBase::GetEdgeNumberOfAdjacentFaces(vtkIdType Edge)
{
	vtkIdList *NonManifoldFaces=this->EdgesNonManifoldFaces[Edge];
	if (NonManifoldFaces!=0)
		return (2+NonManifoldFaces->GetNumberOfIds());

	vtkIdType f1,f2;
	this->GetEdgeFaces(Edge,f1,f2);
	if (f2!=-1)
		return (2);
	if (f1!=-1)
		return (1);
	else
		return (0);
}

vtkIdType vtkSurfaceBase::BisectEdge(vtkIdType Edge)
{
	vtkIdType v1,v2,v3;
	vtkIdList *FList=vtkIdList::New();
	
	this->GetEdgeVertices(Edge,v1,v2);
	this->GetEdgeFaces(Edge,FList);
	
	double P1[3];
	double P2[3];

	this->GetPoint(v1,P1);
	this->GetPoint(v2,P2);
	for (int i=0;i<3;i++)
		P1[i]=0.5*(P1[i]+P2[i]);

	vtkIdType NewVertex=this->AddVertex(P1);
	this->DeleteEdgeInRing(Edge,v1);
	this->Vertex1->SetValue(Edge,NewVertex);
	this->Vertex2->SetValue(Edge,v2);

	for (int i=0;i<FList->GetNumberOfIds();i++)
	{
		vtkIdType Face=FList->GetId(i);
		v3=this->GetThirdPoint(Face,v1,v2);
		this->DeleteFaceInRing(Face,this->IsEdge(v1,v3));
		
		vtkIdType ve1,ve2,ve3;
		this->GetFaceVertices(Face,ve1,ve2,ve3);
		if (((ve1==v1)&&(ve2==v2))
			||((ve2==v1)&&(ve3==v2))
			||((ve3==v1)&&(ve1==v2)))
		{
			this->AddFace(v1,NewVertex,v3);
			this->SetFace(Face,NewVertex,v2,v3);
		}
		else
		{
			this->AddFace(v3,NewVertex,v1);
			this->SetFace(Face,NewVertex,v3,v2);
		}

		this->AddEdge(v3,NewVertex,Face);
	}
	
	this->InsertEdgeInRing(Edge,NewVertex);

	FList->Delete();
	return (NewVertex);
}

vtkIdType FindVertexIndex (vtkIdType *Vertices,vtkIdType Vertex, vtkIdType
NumberOfVertices)
{
	for (vtkIdType i=0;i<NumberOfVertices;i++)
	{
		if (Vertex==Vertices[i])
			return (i);
	}
	return (-1);
}

void vtkSurfaceBase::ChangeFaceVertex(vtkIdType Face, vtkIdType OldVertex, vtkIdType NewVertex)
{
	vtkIdType *Vertices;
	vtkIdType NumberOfVertices;
	this->GetCellPoints(Face,NumberOfVertices,Vertices);

	vtkIdType Index=FindVertexIndex(Vertices,OldVertex,NumberOfVertices);
	assert (Index>=0);

	vtkIdType Neighbours[2];
	Neighbours[0]=Vertices[(Index+1)%NumberOfVertices];
	Neighbours[1]=Vertices[(Index+NumberOfVertices-1)%NumberOfVertices];

	Vertices[Index]=NewVertex;

	for (int i=0;i<2;i++)
	{
		vtkIdType Edge=this->IsEdge(OldVertex,Neighbours[i]);
		if ((Poly2->GetValue(Edge)<0)&&(this->IsEdge(NewVertex,Neighbours[i])<0))
		{
			// this edge has only one adjacent face (the face we are modifying right now).
			// we just modify the edge
			this->Vertex1->SetValue(Edge,NewVertex);
			this->Vertex2->SetValue(Edge,Neighbours[i]);
			this->DeleteEdgeInRing(Edge,OldVertex);
			this->InsertEdgeInRing(Edge,NewVertex);
			this->CleanVertex(OldVertex);
		}
		else
		{
			// this edge has other adjacent faces, we will create a new one.
			this->DeleteFaceInRing(Face,Edge);
			this->CleanEdge(Edge);
			this->AddEdge(NewVertex,Neighbours[i],Face);
		}
	}

	VisitedPolygons->SetValue(Face,0);
	this->ConquerOrientationFromFace(Face);
}

void vtkSurfaceBase::SetOrientationOn()
{
	this->OrientedSurface=true;
}

void vtkSurfaceBase::SetOrientationOff()
{
	this->OrientedSurface=false;
}

bool vtkSurfaceBase::CheckStructure()
{
	bool Problem=false;
	vtkIdType *Vertices,NumberOfVertices;

	for (vtkIdType i=0;i<this->GetNumberOfCells();i++)
	{
		if (this->IsFaceActive(i))
		{
			this->GetFaceVertices(i,NumberOfVertices,Vertices);
			for (vtkIdType j=0;j<NumberOfVertices;j++)
			{
				if (this->IsEdge(Vertices[j],Vertices[(j+1)%NumberOfVertices])<0)
				{
					cout<<"Problem ! Face "<<i<<" misses edge ["<<Vertices[j]<<" "
					<<Vertices[(j+1)%NumberOfVertices]<<"]"<<endl;
					Problem=true;
				}
			}
		}
	}
	return Problem;
}

bool vtkSurfaceBase::IsVertexManifold( const vtkIdType& iV )
{
	vtkIdType v1,v2;
	vtkIdType FirstEdge,FirstVertex;
	vtkIdType f1,f2,f3;
	vtkIdType *Edges,NumberOfRemainingEdges;
	
	this->GetVertexNeighbourEdges(iV,NumberOfRemainingEdges,Edges);

	// if there is only one edge, the vertex is non manifold
	if (NumberOfRemainingEdges<2)
		return (false);
		
	// detect non-manifold edges
	for (vtkIdType i=0;i<NumberOfRemainingEdges;i++)
	{
		if (!this->IsEdgeManifold(Edges[i]))
			return (false);
	}

	FirstEdge=this->GetFirstEdge(iV);
	this->GetEdgeVertices(FirstEdge,FirstVertex,v1);
	if (FirstVertex==iV)
		FirstVertex=v1;
	v1=FirstVertex;

	NumberOfRemainingEdges--;
	this->GetEdgeFaces(FirstEdge,f1,f2);
	
	// turn in the first direction
	v2=this->GetThirdPoint(f1,iV,FirstVertex);
	do
	{
		if (--NumberOfRemainingEdges==0)
			return (true);
		Conquer(f1,iV,v2,f3,v1);
		f1=f3;
		v2=v1;
	} while ((f1>=0)&&(v2!=FirstVertex));

	// if we turned around back to the first vertex, or there is only one direction (f2==-1) return false
	if ((f2<0)||(v2==FirstVertex))
		return (false);

	// turn in the second direction
	v1=FirstVertex;
	v2=this->GetThirdPoint(f2,iV,FirstVertex);
	do
	{
		if (--NumberOfRemainingEdges==0)
			return (true);
		Conquer(f2,iV,v2,f3,v1);
		f2=f3;
		v2=v1;
		
	} while ((f2>=0)&&(v2!=FirstVertex));
	// all adjacent faces have been visited, but there are edges remaining, so return false
	return (false);
}

/// switches the cells orientation (usefull when the mesh is displayed all black...)
void vtkSurfaceBase::SwitchOrientation()
{
	vtkIdType NumberOfVertices, *Vertices;
	this->GetFaceVertices(0,NumberOfVertices,Vertices);
	vtkIdType *Vertices2=new vtkIdType[NumberOfVertices];
	vtkIdType i;
	for (i=0;i<NumberOfVertices;i++)
		Vertices2[i]=Vertices[i];

	for (i=0;i<NumberOfVertices;i++)
		Vertices[i]=Vertices2[NumberOfVertices-1-i];

	this->CheckNormals();	
	delete [] Vertices2;
}

void vtkSurfaceBase::SQueeze()
{
	// Resize edges attributes
	int numEdges=this->GetNumberOfEdges();
	this->Vertex1->Resize(numEdges);
	this->Vertex2->Resize(numEdges);
	this->Poly1->Resize(numEdges);

	int numAllocatedEdges=this->NumberOfAllocatedEdgesAttributes;
	if (numEdges!=numAllocatedEdges)
	{
		if (numEdges>numAllocatedEdges)
			cout<<"Problem while squeezing!!!!!"<<endl;
		vtkIdList **NewArray=0;
		if (numEdges!=0)
		{
			NewArray=new vtkIdList*[numEdges];
			memcpy(NewArray, this->EdgesNonManifoldFaces,numEdges * sizeof(vtkIdList*));
		}
		if (this->EdgesNonManifoldFaces)
		{
			int i;
			for (i=numEdges;i<numAllocatedEdges;i++)
			{
				if (this->EdgesNonManifoldFaces[i])
					this->EdgesNonManifoldFaces[i]->Delete();
			}
			delete [] this->EdgesNonManifoldFaces;
		}
		this->EdgesNonManifoldFaces=NewArray;
	}
	this->NumberOfAllocatedEdgesAttributes=numEdges;

	// Resize vertices Atributes

	vtkIdType numVertices=this->GetNumberOfPoints();
	vtkIdType numAllocatedVertices=this->NumberOfAllocatedVerticesAttributes;

	if (numVertices!=numAllocatedVertices)
	{
		if (numVertices<numAllocatedVertices)
			cout<<"Problem while squeezing!!!!!"<<endl;
		vtkIdType **NewArray=0;
		if (numVertices!=0)
		{
			NewArray=new vtkIdType*[numVertices];
			memcpy(NewArray, this->VerticesAttributes,numVertices* sizeof(int*));
		}
		if (this->VerticesAttributes)
		{
			vtkIdType i;
			for (i=numVertices;i<numAllocatedVertices;i++)
			{
				if (this->VerticesAttributes[i])
					delete[] this->VerticesAttributes[i];
			}
			delete [] this->VerticesAttributes;
		}
		this->VerticesAttributes=NewArray;
	}

	this->NumberOfAllocatedVerticesAttributes=numVertices;

	// Resize polygons Atributes
	this->VisitedPolygons->Resize(this->GetNumberOfCells());

	// Squeeze the PolyData
	this->vtkPolyData::Squeeze();
}

void vtkSurfaceBase::ConquerOrientationFromFace(vtkIdType Face)
{
	std::queue <vtkIdType> EdgesQueue;
	vtkIdType v1,v2,v3,v4,f1,f2,Edge,Edge2;
	vtkIdType Visited1,Visited2;
	vtkIdType NumberOfPoints1,NumberOfPoints2,*Face1,*Face2;
	vtkIdType j;
	vtkIdType FlipFace;

	this->GetCellPoints(Face,NumberOfPoints1,Face1);
	for (j=0;j<NumberOfPoints1;j++)
	{
		Edge=this->IsEdge(Face1[j],Face1[(j+1)%NumberOfPoints1]);
		if (this->IsEdgeManifold(Edge))
			EdgesQueue.push(Edge);
	}

	int Index1,Index2;

	while (EdgesQueue.size())
	{
		// pop an edge from the queue
		Edge=EdgesQueue.front();
		EdgesQueue.pop();
		this->GetEdgeFaces(Edge,f1,f2);
		if (f2>=0)
		{
			Visited1=this->VisitedPolygons->GetValue(f1);
			Visited2=this->VisitedPolygons->GetValue(f2);
			if ((Visited1!=Visited2))
			{
				// if one adjacent face was not visited yet
				if (Visited2==1)
				{
					v1=f1;
					f1=f2;
					f2=v1;
					Visited1=Visited2;
				}
				// f1 was wisited, not f2
				FlipFace=0;
				this->GetCellPoints(f1,NumberOfPoints1,Face1);
				this->GetCellPoints(f2,NumberOfPoints2,Face2);
				this->GetEdgeVertices(Edge,v1,v2);
				Index1=FindVertexIndex(Face1,v1,NumberOfPoints1);
				Index2=FindVertexIndex(Face2,v1,NumberOfPoints2);

				v3=Face1[(Index1+1)%NumberOfPoints1];
				v4=Face2[(Index2+1)%NumberOfPoints2];
				if (((v3==v2)&&(v4==v2))||((v3!=v2)&&(v4!=v2)))
				{
					// flip f2
					for (j=0;j<NumberOfPoints2/2;j++)
					{
						v1=Face2[j];
						Face2[j]=Face2[NumberOfPoints2-1-j];
						Face2[NumberOfPoints2-1-j]=v1;
					}
				}

				this->VisitedPolygons->SetValue(f2,1);
				// add adjacent edges of the face to the queue, except Edge
				for (j=0;j<NumberOfPoints2;j++)
				{
					Edge2=this->IsEdge(Face2[j],Face2[(j+1)%NumberOfPoints2]);
					if ((Edge2!=Edge)&&(this->IsEdgeManifold(Edge2)==1))
						EdgesQueue.push(Edge2);
				}
			}
		}
	}
}

double vtkSurfaceBase::GetValenceEntropy()
{

	int i;		// Loop counter
	int v1;	// Valence of given point
	double inv_number;	// 1.0 / number of points
	double inv_Log2;	// 1.0 / log(2.0)
	double s;			// 
	double p;			//

	// Array of valences, indexed by point's ID
	vtkIntArray *Vals = vtkIntArray::New();
	Vals->Resize(1000);

	// Compute valence for each point and fill Vals
	vtkIdList *List = vtkIdList::New();
	for (i=0; i<1000; i++) Vals->SetValue(i,0);
	for (i=0; i<this->GetNumberOfPoints(); i++)
	{
		this->GetVertexNeighbourFaces(i,List);
		v1 = List->GetNumberOfIds();
		Vals->SetValue(v1,Vals->GetValue(v1)+1);

	}
	List->Delete();

	inv_number = 1.0 / this->GetNumberOfPoints();
	inv_Log2 = 1.0 / log(2.0);
	s = 0;

	for (i=0;i<1000;i++)
	{
		p = Vals->GetValue(i);
		if (p)
		{
			p *= inv_number;
			s -= p * log(p) * inv_Log2;
		}
	}

	Vals->Delete();

	return (s);
}

void vtkSurfaceBase::GetValenceTab(char *filename)
{
	int i;
	int v1;

	int min=0;
	int Max=0;

	// allocate array for valences
	vtkIntArray *Vals=vtkIntArray::New();
	Vals->Resize(1000);
	for (i=0;i<1000;i++) Vals->SetValue(i,0);

	// fill valence's array
	vtkIdList *List=vtkIdList::New();
	for (i=0;i<this->GetNumberOfPoints();i++)
	{
		this->GetVertexNeighbourFaces(i,List);
		v1=List->GetNumberOfIds();
		Vals->SetValue(v1,Vals->GetValue(v1)+1);
	}
	List->Delete();


	// find first non zero value's index
	for(i=0;i<1000;i++)
	{
		if((Vals->GetValue(i)==0)&&(Vals->GetValue(i+1)!=0))
		{
			min = i;
			break;
		}
	}

	// find last non zero value's index
	for(i=0;i<1000;i++)
	{
		if((Vals->GetValue(i)!=0)&&(Vals->GetValue(i+1)==0))
			Max = i+2;
	}

	cout <<min <<endl;
	cout <<Max <<endl;

	std::ofstream out;
	out.open (filename, ofstream::out | ofstream::trunc);

	for(i=min;i<Max;i++)
	{
		if(Vals->GetValue(i)!=0)
			out<<"Degree "<<i<<" : "<<Vals->GetValue(i)+1<<endl;
		else
			out<<"Degree "<<i<<" : 0"<<endl;
	}

	out.close();

	Vals->Delete();

}

void vtkSurfaceBase::GetEdgeFaces(vtkIdType e1,vtkIdList *FList)
{
	FList->Reset();
	vtkIdType Face=this->Poly1->GetValue(e1);
	if (Face==-1)
		return;
	FList->InsertNextId(Face);
	Face=this->Poly2->GetValue(e1);
	if (Face==-1)
		return;
	FList->InsertNextId(Face);
	vtkIdList *FList2=this->EdgesNonManifoldFaces[e1];
	if (FList2==0)
		return;
	int NumberOfFaces=FList2->GetNumberOfIds()-1;
	for (;NumberOfFaces!=-1;NumberOfFaces--)
		FList->InsertNextId(FList2->GetId(NumberOfFaces));
}

/****************************************************************/
/****************************************************************/
void vtkSurfaceBase::DeleteVertex(vtkIdType v1)
{
	this->VerticesGarbage.push(v1);
	this->ActiveVertices->SetValue(v1,0);
}
// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::DeleteEdge(vtkIdType EdgeToRemove)
{
	if (this->Poly1->GetValue(EdgeToRemove)!=-1)
		cout<<"ERROR: Trying to remove an edge which is not free !"<<endl;

	this->ActiveEdges->SetValue(EdgeToRemove,0);
	this->EdgesGarbage.push(EdgeToRemove);
	vtkIdType v1,v2;
	this->GetEdgeVertices(EdgeToRemove,v1,v2);
	this->DeleteEdgeInRing(EdgeToRemove,v1);
	this->DeleteEdgeInRing(EdgeToRemove,v2);
	this->CleanVertex(v1);
	this->CleanVertex(v2);
	this->Vertex2->SetValue(EdgeToRemove,v1);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::DeleteFace(vtkIdType f1)
{
	// test whether the face was already deleted
	if (this->IsFaceActive(f1)==0)
		return;
		
	vtkIdType NumberOfPoints,*Points,Loc;
	Loc = this->Cells->GetCellLocation(f1);
	vtkIdType i;
	vtkIdType Edge;

	this->Polys->GetCell(Loc,NumberOfPoints,Points);
	for (i=0;i<NumberOfPoints;i++)
	{
		Edge=this->IsEdge(Points[i],Points[(i+1)%NumberOfPoints]);
		this->DeleteFaceInRing(f1,Edge);
		this->CleanEdge(Edge);
	}
	for (i=0;i<NumberOfPoints;i++)
	{
		Points[i]=Points[0];
	}

	this->DeleteCell(f1);
	this->CellsGarbage[NumberOfPoints].push(f1);
	this->ActivePolygons->SetValue(f1,0);
}
void vtkSurfaceBase::CleanEdge(vtkIdType Edge)
{
	if (this->CleanEdges==0)
		return;
	if (this->Poly1->GetValue(Edge)==-1)
	{
		this->DeleteEdge(Edge);
	}
}
void vtkSurfaceBase::CleanVertex(vtkIdType Vertex)
{
	if (this->CleanVertices==0)
		return;
	if (this->GetValence(Vertex)==0)
		this->DeleteVertex(Vertex);
}


void vtkSurfaceBase::MergeVertices(vtkIdType v1, vtkIdType v2)
{
	vtkIdList *List=vtkIdList::New();
	vtkIdList *List2=vtkIdList::New();
	
	vtkIdType Edge=this->IsEdge(v1,v2);
	
	// delete polygons adjacent to the edge [v1 v2] (if there are any)
	if (Edge>=0)
	{
		this->GetEdgeFaces(Edge,List);
		for (int i=0;i<List->GetNumberOfIds();i++)
			this->DeleteFace(List->GetId(i));
	}
	
	// for each remaining face adjacent to v2, replace v2 by v1
	this->GetVertexNeighbourFaces(v2,List);
	for (vtkIdType i=0;i<List->GetNumberOfIds();i++)
	{
		vtkIdType NumberOfVertices;
		vtkIdType *Vertices;
		GetFaceVertices(List->GetId(i), NumberOfVertices, Vertices);
		Vertices[FindVertexIndex (Vertices,v2,NumberOfVertices)]=v1;
	}
	
	// Modify every edge adjacent to v2
	this->GetVertexNeighbourEdges(v2,List);
	for (vtkIdType i=0;i<List->GetNumberOfIds();i++)
	{
		vtkIdType v3,v4;
		Edge=List->GetId(i);
		this->GetEdgeVertices(Edge,v3,v4);
		if (v4==v2)
		{
			// ensure that v3==v2
			vtkIdType Vertex=v3;
			v3=v4;
			v4=Vertex;
		}
		
		vtkIdType Edge2;
		Edge2=this->IsEdge(v1,v4);
		if (Edge2>=0)
		{
			// the edge [v1 v4] already exists. Merge it with [v2 v4]
			this->GetEdgeFaces(Edge,List2);
			for (int j=0;j<List2->GetNumberOfIds();j++)
			{
				vtkIdType Face;
				Face=List2->GetId(j);
				this->DeleteFaceInRing(Face,Edge);
				this->InsertFaceInRing(Face,Edge2);
			}
			
			// delete [v2 v4]
			this->DeleteEdge(Edge);
		}
		else
		{
			// the edge [v1 v4] does not exist. Just modify [v2 v4] to [v1 v4]
			this->Vertex1->SetValue(Edge,v1);
			this->Vertex2->SetValue(Edge,v4);

			this->DeleteEdgeInRing(Edge,v2);
			this->InsertEdgeInRing(Edge,v1);
		}
	}
	
	this->DeleteVertex(v2);
	List->Delete();
	List2->Delete();
}

vtkSurfaceBase* vtkSurfaceBase::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkSurfaceBase");
	if(ret)
	{
		return (vtkSurfaceBase*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkSurfaceBase);
}

void vtkSurfaceBase::CheckNormals()
{
	int i;
	for (i=0;i<this->GetNumberOfCells();i++)
		this->VisitedPolygons->SetValue(i,0);

	// We have to loop on all the mesh faces to handle correctly non connnected components
	for (i=0;i<this->GetNumberOfCells();i++)
	{
		if ((this->VisitedPolygons->GetValue(i)==0)&&(this->IsFaceActive(i)==1))
		{
			this->VisitedPolygons->SetValue(i,1);
			this->ConquerOrientationFromFace(i);
		}
	}
	this->OrientedSurface=true;
}

vtkIdType vtkSurfaceBase::IsFace(vtkIdType v1, vtkIdType v2, vtkIdType v3)
{
	vtkIdType edge;
	vtkIdType f1;
	vtkIdType f2;
	int F;
	edge = this->IsEdge(v1,v2);
	
	// test if the edge [v1v2] exists or not
	if (edge<0) 
		return (-1);

	this->GetEdgeFaces(edge,f1,f2);
	
	// test whether the edge [v1v2] is isolated or not
	if (f1<0) 
		return (-1);
		
	if (v3 == this->GetThirdPoint(f1,v1,v2)) 
		return (f1);
	
	// test whether if there is another adjacent face
	if (f2<0) 
		return (-1);

	if (v3==this->GetThirdPoint(f2,v1,v2)) 
		return (f2);

	// test whether there are non-manifold adjacent faces
	vtkIdList *FList2=this->EdgesNonManifoldFaces[edge];
	if (FList2==0)
		return (-1);
	
	int NumberOfFaces=FList2->GetNumberOfIds()-1;
	for (;NumberOfFaces!=-1;NumberOfFaces--)
	{
		F=FList2->GetId(NumberOfFaces);
		if (v3==this->GetThirdPoint(F,v1,v2))
		return (F);
    }
	return (-1);
}

void vtkSurfaceBase :: GetFaceNeighbours(vtkIdType Face,vtkIdListCollection *FList)
{
	vtkIdType v0,v1,v2;
	vtkIdType e0,e1,e2;
		
	vtkIdList *e0Face = vtkIdList :: New();
	vtkIdList *e1Face = vtkIdList :: New();
	vtkIdList *e2Face = vtkIdList :: New();
	
	this->GetFaceVertices(Face,v0,v1,v2);
	e0=this->IsEdge(v0,v1);
	e1=this->IsEdge(v0,v2);
	e2=this->IsEdge(v1,v2);
		
	if(e0!=-1)
	{
		this->GetEdgeFaces(e0,e0Face);
		FList->AddItem(e0Face);
	}

	if(e1!=-1)
	{
		this->GetEdgeFaces(e1,e1Face);
		FList->AddItem(e1Face);
	}

	if(e2!=-1)
	{
		this->GetEdgeFaces(e2,e2Face);
		FList->AddItem(e2Face);
	}
		
	e0Face->Delete();
	e1Face->Delete();
	e2Face->Delete();
}

// ****************************************************************
// ****************************************************************

void vtkSurfaceBase::GetFaceNeighbours(vtkIdType Face,vtkIdList *FList)
{
	vtkIdType NumberOfVertices, *Vertices;
	vtkIdList *List = vtkIdList :: New();
	this->GetFaceVertices(Face, NumberOfVertices, Vertices);
	FList->Reset();

	for (vtkIdType  i=0;i<NumberOfVertices;i++)
	{
		vtkIdType Edge;
		if (i<NumberOfVertices-1)
			Edge=this->IsEdge(Vertices[i],Vertices[i+1]);
		else
			Edge=this->IsEdge(Vertices[i],Vertices[0]);

		this->GetEdgeFaces(Edge,List);
		for (vtkIdType j=0;j<List->GetNumberOfIds();j++)
			FList->InsertUniqueId(List->GetId(j));
	}
	List->Delete();
}

vtkIdType vtkSurfaceBase::FlipEdgeSure(vtkIdType edge)
{
	std::stack<vtkIdType> Edges;
	Edges.push(edge);
	while (!Edges.empty())
	{
		vtkIdType EdgeToFlip=Edges.top();
		vtkIdType ExistingEdge=this->FlipEdge(EdgeToFlip);
		if (ExistingEdge!=-1)
			Edges.push(ExistingEdge);
		else
			Edges.pop();
	}
	return (-1);
}

// ****************************************************************
// ****************************************************************
 vtkIdType vtkSurfaceBase::FlipEdge(vtkIdType edge)
{
	vtkIdType edge1,v1,v2,v3,v4,v7,v8,v9,f1,f2;

	this->GetEdgeFaces(edge,f1,f2);
	if (f2<0)
	{
//		cout<<"*** Problem: flip edge <<"<<edge<<" which is not adjacent to 2 faces !***"<<endl;
		return(edge);
	}
	this->GetEdgeVertices(edge,v1,v2);

	v3=this->GetThirdPoint(f1,v1,v2);
	v4=this->GetThirdPoint(f2,v1,v2);
	edge1=this->IsEdge(v3,v4);
	if (edge1>=0)
	{
//		cout<<"*** Problem: flip edge "<<edge<<" but the resulting edge already exists!***"<<endl;
		return (edge1);
	}

	this->DeleteEdgeInRing(edge,v1);
	this->DeleteEdgeInRing(edge,v2);

	this->Vertex1->SetValue(edge,v3);
	this->Vertex2->SetValue(edge,v4);

	this->InsertEdgeInRing(edge,v3);
	this->InsertEdgeInRing(edge,v4);

	this->GetFaceVertices(f1,v7,v8,v9);
	if (v1==v7)
	{
		if (v2==v8)
		{
			this->SetFace(f1,v7,v4,v9);
			this->SetFace(f2,v4,v8,v9);
		}
		else
		{
			this->SetFace(f1,v7,v8,v4);
			this->SetFace(f2,v4,v8,v9);
		}
	}
	else
	{
		if (v1==v8)
		{
			if (v2==v9)
			{
				this->SetFace(f1,v8,v4,v7);
				this->SetFace(f2,v7,v4,v9);
			}
			else
			{
				this->SetFace(f1,v9,v4,v8);
				this->SetFace(f2,v4,v9,v7);
			}
		}
		else
		{
			if (v2==v7)
			{
				this->SetFace(f1,v1,v4,v8);
				this->SetFace(f2,v8,v4,v7);
			}
			else
			{
				this->SetFace(f1,v9,v3,v4);
				this->SetFace(f2,v8,v4,v3);
			}
		}
	}

	edge1=this->IsEdge(v1,v4);
	if (this->Poly1->GetValue(edge1)==f2)
		this->Poly1->SetValue(edge1,f1);
	else
		this->Poly2->SetValue(edge1,f1);

	edge1=this->IsEdge(v2,v3);
	if (this->Poly1->GetValue(edge1)==f1)
		this->Poly1->SetValue(edge1,f2);
	else
		this->Poly2->SetValue(edge1,f2);

	return (-1);
}

// ** METHODE GetNumberOfBoundaries
int vtkSurfaceBase::GetNumberOfBoundaries(vtkIdType v1)
{
	vtkIdType i;

	int numBoundaries = 0;
	vtkIdType *Edges,NumberOfEdges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);

	for (i=0;i<NumberOfEdges;i++)
	{
		if (Poly2->GetValue(Edges[i])<0)
			numBoundaries++;
	}
	if (numBoundaries % 2 == 1)
	{
		//printf("numBoundaries = %d... what?\n", numBoundaries);
		return (numBoundaries + 1) / 2;
	}
	else
		return (numBoundaries/2);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetVertexNeighbourEdges(vtkIdType v1, vtkIdList *Output)
{
	Output->Reset();
	vtkIdType *Edges,NumberOfEdges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	for (vtkIdType i=0;i<NumberOfEdges;i++)
		Output->InsertNextId(Edges[i]);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetVertexNeighbourFaces(vtkIdType v1, vtkIdList *Output)
{
	vtkIdType f1;
	vtkIdType f2;
	vtkIdType NumberOfEdges,*Edges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	Output->Reset();
	int i;
	for (i=0;i<NumberOfEdges;i++)
	{
		this->GetEdgeFaces(Edges[i],f1,f2);
		if (f1>=0) Output->InsertUniqueId(f1);
		if (f2>=0) Output->InsertUniqueId(f2);
		
		vtkIdList *OtherFaces=this->EdgesNonManifoldFaces[Edges[i]];
		if (OtherFaces)
		{
			for (int j=0;j<OtherFaces->GetNumberOfIds();j++)
				Output->InsertUniqueId(OtherFaces->GetId(j));
		}
	}
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetVertexNeighbours(vtkIdType v1, vtkIdList *Output)
{
	vtkIdType v2,edge;
	vtkIdType i;
	Output->Reset();
	vtkIdType NumberOfEdges,*Edges;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	for (i=0;i<NumberOfEdges;i++)
	{
		edge=Edges[i];

		v2=this->Vertex1->GetValue(edge);

		if (v2==v1)	Output->InsertNextId(this->Vertex2->GetValue(edge));
		else		Output->InsertNextId(v2);
	}
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::GetNeighbours(vtkIdList *Input,vtkIdList *Output)
{
	vtkIdType i,v1,v2,edge,j;
	Output->Reset();
	vtkIdType NumberOfEdges,*Edges;

	for (i=0;i<Input->GetNumberOfIds();i++)
	{
		v1 = Input->GetId(i);
		Output->InsertUniqueId(v1);
		this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
		for (j=0;j<NumberOfEdges;j++)
		{
			edge = Edges[j];

			v2 = this->Vertex1->GetValue(edge); 

			if (v2==v1) Output->InsertUniqueId(this->Vertex2->GetValue(edge));
			else		Output->InsertUniqueId(v2);
		}
	}
}


vtkIdType vtkSurfaceBase::IsEdgeBetweenFaces(vtkIdType f1, vtkIdType f2)
{
	vtkIdType *Vertices1;
	vtkIdType NumberOfVertices1;
	this->GetCellPoints(f1,NumberOfVertices1,Vertices1);
	vtkIdType *Vertices2;
	vtkIdType NumberOfVertices2;
	this->GetCellPoints(f2,NumberOfVertices2,Vertices2);
	
	vtkIdType Vertices[2];
	
	int Count=0;
	for (vtkIdType i=0;i<NumberOfVertices1;i++)
	{
		vtkIdType v1=Vertices1[i];
		for (vtkIdType j=0;j<NumberOfVertices2;j++)
		{
			vtkIdType v2=Vertices2[j];
			if (v1==v2)
			{
				Vertices[Count++]=v1;
				break;
			}
		}
		if (Count==2)
			return (this->IsEdge(Vertices[0],Vertices[1]));		
	}
	return (-1);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::InsertEdgeInRing(vtkIdType e1,vtkIdType v1)
{
	vtkIdType NumberOfEdges,*Edges;
	vtkIdType *VertexAttributes=this->VerticesAttributes[v1];
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	
	// test wether the vertices were allocated
	if (v1>=this->NumberOfAllocatedVerticesAttributes)
		this->AllocateVerticesAttributes(v1+1);
	vtkIdType  *NewArray;

	vtkIdType  i;
	// First test if the edge is not already in ring
	for (i=0;i<NumberOfEdges;i++)
	{
		if (Edges[i]==e1)
			return;
	}

	// allocate new ring if needed
	if (NumberOfEdges>=VertexAttributes[VERTEX_NUMBER_OF_EDGES_SLOTS])
	{
		NewArray=new vtkIdType[NumberOfEdges+1+VERTEX_EDGES];
		for (i=0;i<NumberOfEdges+VERTEX_EDGES;i++)
		{
			NewArray[i]=VertexAttributes[i];
		}
		NewArray[VERTEX_NUMBER_OF_EDGES_SLOTS]=NumberOfEdges+1;

		delete [] VertexAttributes;
		Edges=NewArray+VERTEX_EDGES;
		this->VerticesAttributes[v1]=NewArray;
		VertexAttributes=NewArray;
	}
	Edges[NumberOfEdges]=e1;
	VertexAttributes[VERTEX_NUMBER_OF_EDGES]=NumberOfEdges+1;
}

void vtkSurfaceBase::DeleteEdgeInRing(vtkIdType e1,vtkIdType v1)
{
	vtkIdType  NumberOfEdges,*Edges;
	vtkIdType  i;
	this->GetVertexNeighbourEdges(v1,NumberOfEdges,Edges);
	for (i=0;i<NumberOfEdges;i++)
	{
		if (Edges[i]==e1)
			break;
	}
	if (i<NumberOfEdges-1)
		Edges[i]=Edges[NumberOfEdges-1];
	this->VerticesAttributes[v1][VERTEX_NUMBER_OF_EDGES]--;
}

void vtkSurfaceBase::AllocateMoreVerticesAttributes()
{
	int NumberOfAttributes=this->NumberOfAllocatedVerticesAttributes;
	double number;
	int NewNumberOfAttributes;
	number=NumberOfAttributes;
	number=1.0+number*1.1;
	NewNumberOfAttributes=(int) number;
	this->AllocateVerticesAttributes(NewNumberOfAttributes);
}

void vtkSurfaceBase::AllocateMoreEdgesAttributes()
{
	int NumberOfAttributes=this->NumberOfAllocatedEdgesAttributes;
	double number;
	int NewNumberOfAttributes;
	number=NumberOfAttributes;
	number=1.0+number*1.1;
	NewNumberOfAttributes=(int) number;
	this->AllocateEdgesAttributes(NewNumberOfAttributes);
}
void vtkSurfaceBase::AllocateMorePolygonsAttributes()
{
	int NumberOfAttributes=this->NumberOfAllocatedPolygonsAttributes;
	double number;
	int NewNumberOfAttributes;
	number=NumberOfAttributes;
	number=1.0+number*1.1;
	NewNumberOfAttributes=(int)number;
	this->AllocatePolygonsAttributes(NewNumberOfAttributes);
}

// ****************************************************************
// ****************************************************************
void vtkSurfaceBase::InsertFaceInRing(vtkIdType Face,vtkIdType Edge)
{
	vtkIdType F1,F2;
	vtkIdList *FList;
	this->GetEdgeFaces(Edge,F1,F2);
	if (F1==-1)
	{
		this->Poly1->SetValue(Edge,Face);
		this->Poly2->SetValue(Edge,-1);
		return;
	}
	if (F1==Face)
		return;

	if (F2==-1)
	{
		this->Poly2->SetValue(Edge,Face);
		return;
	}
	if (F2==Face)
		return;

	FList=this->EdgesNonManifoldFaces[Edge];
	if (!FList)
	{
		FList=vtkIdList::New();
		this->EdgesNonManifoldFaces[Edge]=FList;
	}
	FList->InsertUniqueId(Face);
	return;
}

void vtkSurfaceBase::DeleteFaceInRing(vtkIdType Face,vtkIdType Edge)
{
	vtkIdList *FList=this->EdgesNonManifoldFaces[Edge];
	vtkIdType *F1=this->Poly1->GetPointer(Edge);
	vtkIdType *F2=this->Poly2->GetPointer(Edge);
	if (FList==0)
	{
		if (*F1==Face)
			*F1=*F2;
		*F2=-1;
		return;
	}
	vtkIdType Index=FList->IsId(Face);
	if (Index>=0)
		FList->DeleteId(Face);        
	else
	{
		int LastPosition=FList->GetNumberOfIds()-1;
		vtkIdType F3=FList->GetId(LastPosition);
		FList->DeleteId(F3);
		if (*F1==Face)
			*F1=F3;
		else
			*F2=F3;
	}
	if (FList->GetNumberOfIds()==0)
	{
		FList->Delete();
		this->EdgesNonManifoldFaces[Edge]=0;
	}
}

vtkIdType vtkSurfaceBase::AddVertex(double x, double y, double z)
{
	vtkIdType v1;
	if (this->VerticesGarbage.empty())
	{
		v1=this->GetPoints()->InsertNextPoint(x,y,z);
		if (v1>=this->NumberOfAllocatedVerticesAttributes)
		{
			this->AllocateMoreVerticesAttributes();
		}
		this->VerticesAttributes[v1][VERTEX_NUMBER_OF_EDGES]=0;
		this->ActiveVertices->SetValue(v1,1);
		return (v1);
	}
	else
	{
		v1=this->VerticesGarbage.front();
		this->VerticesGarbage.pop();
		this->VerticesAttributes[v1][VERTEX_NUMBER_OF_EDGES]=0;
		this->GetPoints()->SetPoint(v1,x,y,z);
		this->ActiveVertices->SetValue(v1,1);
		return (v1);
	}
}

// ****************************************************************
// ****************************************************************
vtkIdType vtkSurfaceBase::AddEdge(vtkIdType v1,vtkIdType v2,vtkIdType f1)
{
	while (v1>=this->NumberOfAllocatedVerticesAttributes)
		this->AllocateMoreVerticesAttributes();
	while (v2>=this->NumberOfAllocatedVerticesAttributes)
		this->AllocateMoreVerticesAttributes();
	vtkIdType edge=this->IsEdge(v1,v2);
	if (v1==v2)
	{
		cout<<"Error : creation of a self-loop for vertex "<<v1<<endl;
		return (-1);
	}

	if (edge>=0)
	{
		if (Poly1->GetValue(edge)<0)
		{
			Poly1->SetValue(edge,f1);
			return (edge);
		}

		if (Poly2->GetValue(edge)>=0)
		{
			vtkIdList *LIST=this->EdgesNonManifoldFaces[edge];
			if (!LIST)
			{
				LIST=vtkIdList::New();
				this->EdgesNonManifoldFaces[edge]=LIST;
			}
			LIST->InsertNextId(f1);
			return (edge);
		}
		else
		{
			Poly2->SetValue(edge,f1);
			return (edge);
		}
	}

	if (this->EdgesGarbage.empty())
	{
		if (this->NumberOfEdges>=this->NumberOfAllocatedEdgesAttributes)
			this->AllocateMoreEdgesAttributes();

		edge=this->NumberOfEdges++;
	}
	else
	{
		edge=this->EdgesGarbage.front();
		this->EdgesGarbage.pop();
	}
	this->Vertex1->SetValue(edge,v1);//
	this->Vertex2->SetValue(edge,v2);//
	this->Poly1->SetValue(edge,f1);
	this->Poly2->SetValue(edge,-1);
	this->EdgesNonManifoldFaces[edge]=0;
	this->ActiveEdges->SetValue(edge,1);

	this->InsertEdgeInRing(edge,v1);
	this->InsertEdgeInRing(edge,v2);
	return (edge);
}

vtkIdType vtkSurfaceBase::AddFace(const vtkIdType& v1,const vtkIdType& v2,const
vtkIdType& v3)
{
	int Number=3;
	vtkIdType Vertices[3];
	Vertices[0]=v1;
	Vertices[1]=v2;
	Vertices[2]=v3;
	return this->AddPolygon(Number,Vertices);
}

// ****************************************************************
// ****************************************************************
vtkIdType vtkSurfaceBase::AddPolygon(int NumberOfVertices,vtkIdType *Vertices)
{
	vtkIdType face;
	int Loc;
	int i;

	int Type=0;
	
	// Check if the vertices of the polygon have been created
	
	for (i=0;i<NumberOfVertices;i++)
	{
		if (Vertices[i]>=this->GetNumberOfPoints())
			Type=99;		
	}
	if (Type)
	{
		cout<<"ERROR : attempt to create the Cell with vertices : ";
		for (i=0;i<NumberOfVertices;i++)
			cout<<Vertices[i]<<" ";
		cout<<endl;
		cout<<"But Only "<<this->GetNumberOfPoints()<<" vertices have been created. Exiting program"<<endl;
		exit(1);
	}

	
	switch (NumberOfVertices)
	{
	case 3:
		Type=VTK_TRIANGLE;
		break;
	case 4:
		Type=VTK_QUAD;
		break;

	default:
		Type=VTK_POLYGON;
		break;
	}

	if (this->CellsGarbage[NumberOfVertices].empty())
	{
		face=this->Polys->InsertNextCell(NumberOfVertices);
		for (i=0;i<NumberOfVertices;i++)
			this->Polys->InsertCellPoint(Vertices[i]);
		Loc=this->Polys->GetInsertLocation(NumberOfVertices);
		this->Cells->InsertNextCell(Type,Loc);
		while (face>=this->NumberOfAllocatedPolygonsAttributes)
			this->AllocateMorePolygonsAttributes();
	}
	else
	{
		face=this->CellsGarbage[NumberOfVertices].front();
		this->CellsGarbage[NumberOfVertices].pop();
		Loc = this->Cells->GetCellLocation(face);
		this->Cells->InsertCell(face,Type,Loc);
		this->Polys->ReplaceCell(Loc,NumberOfVertices,Vertices);
	}

	for (i=0;i<NumberOfVertices;i++)
	{
		this->AddEdge(Vertices[i],Vertices[(i+1)%NumberOfVertices],face);
	}

	if ((this->FirstTime)&&(face==0))
	{
		this->VisitedPolygons->SetValue(face,1);
		this->FirstTime=false;
	}
	else
		this->VisitedPolygons->SetValue(face,0);

	this->ActivePolygons->SetValue(face,1);
	if (this->OrientedSurface)
	{
		this->ConquerOrientationFromFace(face);
	}
	return (face);
}

// ****************************************************************
// ****************************************************************


void vtkSurfaceBase::AllocateVerticesAttributes(int NumberOfVertices)
{
	if (this->NumberOfAllocatedVerticesAttributes>=NumberOfVertices)
		return;

	vtkIdType **NewArray=new vtkIdType*[NumberOfVertices];
	vtkIdType i;

	for (i=0;i<this->NumberOfAllocatedVerticesAttributes;i++)
		NewArray[i]=this->VerticesAttributes[i];
	for (i=this->NumberOfAllocatedVerticesAttributes;i<NumberOfVertices;i++)
	{
		NewArray[i]=new vtkIdType[VERTEX_EDGES+6];
		NewArray[i][VERTEX_NUMBER_OF_EDGES]=0;
		NewArray[i][VERTEX_NUMBER_OF_EDGES_SLOTS]=6;
	}
	if (this->VerticesAttributes)
		delete []this->VerticesAttributes;
	this->VerticesAttributes=NewArray;
	this->NumberOfAllocatedVerticesAttributes=NumberOfVertices;

	if (!this->ActiveVertices)
		this->ActiveVertices=vtkBitArray::New();
	this->ActiveVertices->Resize(NumberOfVertices);

}

void vtkSurfaceBase::AllocatePolygonsAttributes(int NumberOfPolygons)
{
	if (this->NumberOfAllocatedPolygonsAttributes>=NumberOfPolygons)
		return;
	if (!this->VisitedPolygons)
		this->VisitedPolygons=vtkBitArray::New();
	this->VisitedPolygons->Resize(NumberOfPolygons);
	this->NumberOfAllocatedPolygonsAttributes=NumberOfPolygons;

	if (!this->ActivePolygons)
		this->ActivePolygons=vtkBitArray::New();
	this->ActivePolygons->Resize(NumberOfPolygons);

}

void vtkSurfaceBase::AllocateEdgesAttributes(int NumberOfEdges)
{
	if (this->NumberOfAllocatedEdgesAttributes>=NumberOfEdges)
		return;

	// create a new tab for non manifold edges
	int EDGE;
	vtkIdList **NewArray=new vtkIdList*[NumberOfEdges];
	// copy the old tab into the new one
	for (EDGE=0;EDGE<this->NumberOfEdges;EDGE++)
		NewArray[EDGE]=this->EdgesNonManifoldFaces[EDGE];
	for (EDGE=this->NumberOfEdges;EDGE<NumberOfEdges;EDGE++)
		NewArray[EDGE]=0;
	//replace the old tab by the new
	if (this->EdgesNonManifoldFaces)
		delete [] this->EdgesNonManifoldFaces;
	this->EdgesNonManifoldFaces=NewArray;

	if (!this->Poly2)
		this->Poly1=vtkIdTypeArray::New();
	if (!this->Poly2)
		this->Poly2=vtkIdTypeArray::New();
	if (!this->Vertex1)
		this->Vertex1=vtkIdTypeArray::New();
	if (!this->Vertex2)
		this->Vertex2=vtkIdTypeArray::New();

	this->Poly1->Resize(NumberOfEdges);
	this->Poly2->Resize(NumberOfEdges);
	this->Vertex1->Resize(NumberOfEdges);
	this->Vertex2->Resize(NumberOfEdges);

	if (!this->ActiveEdges)
		this->ActiveEdges=vtkBitArray::New();
	this->ActiveEdges->Resize(NumberOfEdges);

	this->NumberOfAllocatedEdgesAttributes=NumberOfEdges;
}
// ****************************************************************
// ****************************************************************
// cette fonction initialise l'objet a partir d'un autre object de meme type
// Alloue la place et l'objet pour stocker le meme nombre de points de
// triangles et d'arretes , MAIS NE COPIE PAS LES DONNEES
void vtkSurfaceBase::Init(vtkSurfaceBase *mesh)
{

	// recupere le nombre de triangles, points et arretes
	int numPoints=mesh->GetNumberOfPoints();
	int numFaces=mesh->GetNumberOfCells();
	int numEdges=mesh->GetNumberOfEdges();

	// delegue le travail
	this->Init(numPoints,numFaces,numEdges);

}


// ****************************************************************
// ****************************************************************
// cette fonction initialise l'objet a partir d'un nombre de points,
// de triangles et d'arretes , renvoie un objet vide
void vtkSurfaceBase::Init(int numPoints, int numFaces, int numEdges)
{
	vtkCellArray *CellsArray1;
	CellsArray1 = vtkCellArray::New();
	CellsArray1->Allocate(4*numFaces,numFaces);
	this->SetPolys(CellsArray1);
	CellsArray1->Delete();

	if (this->Cells)
		this->Cells->Delete();
	
	this->Cells = vtkCellTypes::New();
	this->Cells->Allocate(numFaces,numFaces);

	// create and allocate memory for points
	vtkPoints *Points1=vtkPoints::New();
	Points1->Allocate(numPoints);
	this->SetPoints(Points1);
	Points1->Delete();

	// create and allocate memory for all the vtkSurfaceBase specific tables
	this->AllocateVerticesAttributes(numPoints);
	this->AllocateEdgesAttributes(numEdges);
	this->AllocatePolygonsAttributes(numFaces);
}

// ****************************************************************
// ****************************************************************
// fonction CreateFromPolyData
// cette fonction initialise l'objet a partir d'un polydata
// copie les donnees (points et triangles)
// et cree le tableau d'arretes
void vtkSurfaceBase::CreateFromPolyData(vtkPolyData *input)
{
	vtkIdType i,j,v1,v2;
	vtkIdType NumberOfVertices,*Vertices;

	// just copy the polydata in input
	this->ShallowCopy(input);

	// Delete the cells that are not polygons
	vtkCellArray *VerticesCells=this->GetVerts();
	vtkCellArray *LinesCells=this->GetLines();
	vtkCellArray *StripsCells=this->GetStrips();

	if ((VerticesCells!=0)||(LinesCells!=0)||(StripsCells!=0))
	{
		this->SetVerts(0);
		this->SetLines(0);
		this->SetStrips(0);
		this->Modified();
		this->BuildCells();
	}

	vtkIdType numPoints=this->GetNumberOfPoints();
	vtkIdType numFaces=this->GetNumberOfCells();

	this->AllocateVerticesAttributes(numPoints);
	this->AllocateEdgesAttributes(numPoints+numFaces+1000);
	this->AllocatePolygonsAttributes(numFaces);

	for (i=0;i<this->GetNumberOfCells();i++)
	{
		bool ActiveFace=false;
		this->GetCellPoints(i,NumberOfVertices,Vertices);
		if (NumberOfVertices>1)
		{
			// test whether the face is an active one (at least its first two vertices should be different)
			if (Vertices[0]!=Vertices[1])
			{
				ActiveFace=true;
				for (j=0;j<NumberOfVertices;j++)
				{	
					v1=Vertices[j];
					v2=Vertices[(j+1)%NumberOfVertices];
					this->AddEdge(v1,v2,i);
				}
			}
		}
		if (ActiveFace)
			this->ActivePolygons->SetValue(i,1);
		else
		{
			// The face is not used. Let's push it in the garbage collector
			this->ActivePolygons->SetValue(i,0);
			this->CellsGarbage[NumberOfVertices].push(i);			
		}
	}

	if (this->OrientedSurface)
	{
		this->CheckNormals();
	}
}

vtkSurfaceBase::vtkSurfaceBase()
{
	this->FirstTime=true;
	this->SetOrientationOn();

	this->NumberOfEdges=0;
	this->NumberOfAllocatedVerticesAttributes=0;
	this->NumberOfAllocatedEdgesAttributes=0;
	this->NumberOfAllocatedPolygonsAttributes=0;

	// vertices attributes
	this->VerticesAttributes=0;
	this->ActiveVertices=0;


	// edges attributes
	this->Poly1 = 0;
	this->Poly2 = 0;
	this->Vertex1 = 0;
	this->Vertex2 = 0;
	this->EdgesNonManifoldFaces=0;
	this->ActiveEdges=0;

	// faces attributes 
	this->VisitedPolygons=0;
	this->ActivePolygons=0;

	this->CleanEdges=1;
	this->CleanVertices=0;
	
	this->Init(50,100,150);
}

// ****************************************************************
// ****************************************************************
vtkSurfaceBase::~vtkSurfaceBase() //Destructeur
{

	vtkIdType i;
	if (this->VerticesAttributes)
	{
		vtkIdType *Ring;
		for (i=0;i<this->NumberOfAllocatedVerticesAttributes;i++)
		{
			Ring=this->VerticesAttributes[i];
			if (Ring)
				delete [] Ring;
		}
		delete [] this->VerticesAttributes;
	}

	if (this->Poly1) this->Poly1->Delete();
	if (this->Poly2) this->Poly2->Delete();
	if (this->Vertex1) this->Vertex1->Delete();
	if (this->Vertex2) this->Vertex2->Delete();

	if (this->EdgesNonManifoldFaces)
	{
		vtkIdList*List;
		for (i=0;i<this->NumberOfAllocatedEdgesAttributes;i++)
		{
			List=this->EdgesNonManifoldFaces[i];
			if (List)
				List->Delete();
		}
		delete [] this->EdgesNonManifoldFaces;
	}
	if (this->VisitedPolygons)
		this->VisitedPolygons->Delete();
		
	if (this->ActivePolygons)
		this->ActivePolygons->Delete();
		
	if (this->ActiveEdges)
		this->ActiveEdges->Delete();
		
	if (this->ActiveVertices)
		this->ActiveVertices->Delete();
}
