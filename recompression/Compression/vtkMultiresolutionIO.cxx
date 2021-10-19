/*=========================================================================
	Program:   Mailleur 3D multi-resolution (Creatis 1997 ~)
	Module:    vtkMultiresolutionIO.cxx
	Language:  C++
	Date:      2003/05
	Auteurs:   Sebastien Valette
	This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME vtkMultiresolutionIO
// .SECTION Description

#include <list>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkLookupTable.h>
#include <vtkWindowToImageFilter.h>
#include <vtkBMPWriter.h>
#include <vtkPLYWriter.h>
#include <vtkCamera.h>

#include "vtkSurface.h"
#include "vtkMultiresolutionIO.h"

// for ROI
#include "vtkSurfaceIterators.h"

// for obj
#include <vtkOBJExporter.h>

// for random choice
#include <iostream>
#include <random>

#include <time.h>
#define MIN(a,b) (((a)<(b))?(a):(b))


void vtkMultiresolutionIO::PrintInsignificantCoeffs()
{
	int coord=0;
	std::list<WaveletCoefficient>::iterator CoeffIter;

	cout<<"InsignificantCoeffs:"<<endl;
	for (CoeffIter=this->InsignificantCoeffs[coord].begin();
		CoeffIter!=this->InsignificantCoeffs[coord].end();
		CoeffIter++)
	{
		cout<<"Filter "<<CoeffIter->FilterId<<", Wavelet "<<CoeffIter->WaveletId<<endl;
	}
}
void vtkMultiresolutionIO::PrintWavelets()
{
	int i,j,count=1;
	for	(i=0;i<this->NumberOfFilters;i++)
	{
		for (j=0;j<this->Filters[i]->Wavelets->GetNumberOfTuples();j++)
		{
			double	value1,value2,value3,Wav[3];
			this->Filters[i]->Wavelets->GetTuple(j,Wav);
			value1=Wav[0];
			value2=Wav[1];
			value3=Wav[2];
			cout<<count++<<" : Filter "<<i<<", wavelet "<<j<<" : "<<value1<<" "<<value2<<" "<<value3<<endl;
		}
	}
}
void vtkMultiresolutionIO::SaveWavelets()
{
	int i,j;
	fstream	WavF;
	WavF.open ("wavelets.dat", ofstream::out | ofstream::trunc|ios::binary);
	for	(i=0;i<this->NumberOfFilters;i++)
	{
		for (j=0;j<this->Filters[i]->Wavelets->GetNumberOfTuples();j++)
		{
			double	value1,value2,value3,Wav[3];
			this->Filters[i]->Wavelets->GetTuple(j,Wav);
			value1=Wav[0];
			value2=Wav[1];
			value3=Wav[2];
			WavF.write((char*) &value1,sizeof(double));
			WavF.write((char*) &value2,sizeof(double));
			WavF.write((char*) &value3,sizeof(double));
		}
	}
	WavF.close();
}

void vtkMultiresolutionIO::LoadWavelets()
{
	fstream WavF;
	WavF.open ("wavelets.dat", ofstream::in|ios::binary);
	int i,j;

	for	(i=0;i<this->NumberOfFilters;i++)
	{
		for (j=0;j<this->Filters[i]->Wavelets->GetNumberOfTuples();j++)
		{
			double	value1,value2,value3,Wav[3];
			WavF.read((char*) &value1,sizeof	(double));
			WavF.read((char*) &value2,sizeof	(double));
			WavF.read((char*) &value3,sizeof	(double));
			Wav[0]=value1;
			Wav[1]=value2;
			Wav[2]=value3;
			this->Filters[i]->Wavelets->SetTuple(j,Wav);
		}
	}
	WavF.close();
}

void vtkMultiresolutionIO::EncodeLiftingProperties()
{
	if (this->Lifting==0)
		this->ArithmeticCoder->EncodeByte(1);
	else
	{
		if (this->Lifting==2)
			this->ArithmeticCoder->EncodeByte(0);
		else
			this->ArithmeticCoder->EncodeByte(this->LiftingRadius+2);
	}
	this->ArithmeticCoder->EncodeBit(this->GeometryPrediction!=0);
}
void vtkMultiresolutionIO::DecodeLiftingProperties()
{
	int LiftingTest=this->ArithmeticCoder->DecodeByte();
	LiftingTest-=2;
	if (LiftingTest>=0)
	{
		this->Lifting=1;
		this->LiftingRadius=LiftingTest;
	}
	else
	{
		if (LiftingTest==-1)
		{
			this->Lifting=0;
			this->LiftingRadius=0;			
		}
		else
		{
			this->Lifting=2;
			this->LiftingRadius=0;
		}
	}
	this->GeometryPrediction=(this->ArithmeticCoder->DecodeBit()==true);

}

void vtkMultiresolutionIO::EncodeMeshConnectivity(vtkSurface *Mesh)
{
	vtkIdType i;
	vtkIdType n,v1,v2,v3;

	this->ArithmeticCoder->EncodeInt(Mesh->GetNumberOfPoints());	// Me: Word -> Int
	this->ArithmeticCoder->EncodeInt(Mesh->GetNumberOfCells());	// Me: Word -> Int
	qsmodel QsmConnectivity;
	QsmConnectivity.initqsmodel(Mesh->GetNumberOfPoints(),23,1000,NULL,1);

	n=Mesh->GetNumberOfCells();
	for (i=0;i<n;i++)
	{
		Mesh->GetFaceVertices(i,v1,v2,v3);
		this->ArithmeticCoder->Encode(v1,&QsmConnectivity);
		this->ArithmeticCoder->Encode(v2,&QsmConnectivity);
		this->ArithmeticCoder->Encode(v3,&QsmConnectivity);
	}
}
vtkSurface *vtkMultiresolutionIO::DecodeMeshConnectivity()
{
	vtkIdType i;
	int v1,v2,v3;
	int NumberOfPoints,NumberOfFaces;
	double P[3];
	P[0]=0;
	P[1]=0;
	P[2]=0;

	vtkSurface *Mesh=vtkSurface::New();
	NumberOfPoints=this->ArithmeticCoder->DecodeInt();	// Me: Word -> Int
	NumberOfFaces=this->ArithmeticCoder->DecodeInt();	// Me: Word -> Int
	Mesh->Init(NumberOfPoints,NumberOfFaces,NumberOfPoints+NumberOfFaces+1000);


	for (i=0;i<NumberOfPoints;i++)
		Mesh->AddVertex(P);

	qsmodel QsmConnectivity;
	QsmConnectivity.initqsmodel(NumberOfPoints,23,1000,NULL,0);


	for (i=0;i<NumberOfFaces;i++)
	{
		v1=this->ArithmeticCoder->Decode(&QsmConnectivity);
		v2=this->ArithmeticCoder->Decode(&QsmConnectivity);
		v3=this->ArithmeticCoder->Decode(&QsmConnectivity);
		Mesh->AddFace(v1,v2,v3);
	}

	vtkSurface *Mesh2=vtkSurface::New();
	Mesh2->CreateFromPolyData(Mesh);
	Mesh->Delete();
	return (Mesh2);


}

void vtkMultiresolutionIO::EncodeMeshGeometry(vtkSurface *Mesh)
{

	vtkIdType i,j;
	double P[3];
	int NumberOfPoints;

	NumberOfPoints=Mesh->GetNumberOfPoints();

	if (this->ArithmeticType==0)
	{
		for (i=0;i<NumberOfPoints;i++)
		{
			Mesh->GetPointCoordinates(i,P);
			for (j=0;j<3;j++)
			{
				this->ArithmeticCoder->EncodeFloat(P[j]);
			}
		}
	}
	else
	{
		int Value;
		for (i=0;i<NumberOfPoints;i++)
		{
			Mesh->GetPointCoordinates(i,P);
			for (j=0;j<3;j++)
			{
				Value=(int) P[j];
				this->ArithmeticCoder->EncodeInt(Value+ 1073741824);	// Me: Word -> Int
			}
		}
	}
}
void vtkMultiresolutionIO::DecodeMeshGeometry(vtkSurface *Mesh)
{

	vtkIdType i,j;
	double P[3];
	int NumberOfPoints;

	NumberOfPoints=Mesh->GetNumberOfPoints();
	if (this->ArithmeticType==0)
	{
		for (i=0;i<NumberOfPoints;i++)
		{
			for (j=0;j<3;j++)
			{		
				P[j]=this->ArithmeticCoder->DecodeFloat();
			}
			Mesh->SetPointCoordinates(i,P);
		}
	}
	else
	{
		for (i=0;i<NumberOfPoints;i++)
		{
			for (j=0;j<3;j++)
			{		
				P[j]=(this->ArithmeticCoder->DecodeInt()- 1073741824);	// Me: Word -> Int
			}
			Mesh->SetPointCoordinates(i,P);
		}
	}
}

void vtkMultiresolutionIO::EncodeScalingFactors(vtkSurface *Mesh)
{
	double Factor, Tx, Ty, Tz;
	Mesh->GetScalingFactors(Factor,Tx,Ty,Tz);
	this->ArithmeticCoder->EncodeFloat((float) Factor);
	this->ArithmeticCoder->EncodeFloat((float) Tx);
	this->ArithmeticCoder->EncodeFloat((float) Ty);
	this->ArithmeticCoder->EncodeFloat((float) Tz);
}
void vtkMultiresolutionIO::DecodeScalingFactors(vtkSurface *Mesh)
{
	double Factor, Tx, Ty, Tz;
	Factor=this->ArithmeticCoder->DecodeFloat();
	Tx=this->ArithmeticCoder->DecodeFloat();
	Ty=this->ArithmeticCoder->DecodeFloat();
	Tz=this->ArithmeticCoder->DecodeFloat();
	Mesh->SetScalingFactors(Factor,Tx,Ty,Tz);
}


void vtkMultiresolutionIO::SetLifting(int lifting)
{
	int i;
	this->Lifting=lifting;

	if (this->NumberOfFilters>0)
		for (i=0;i<this->NumberOfFilters;i++)
			this->Filters[i]->SetLifting(lifting);
}
void vtkMultiresolutionIO::SetArithmeticType (int Type)
{
	int i;
	this->ArithmeticType=Type;
	if (this->NumberOfFilters>0)
		for (i=0;i<this->NumberOfFilters;i++)
			this->Filters[i]->SetArithmeticType(Type);
};

void vtkMultiresolutionIO::SetLiftingRadius(int radius)
{
	int i;
	this->LiftingRadius=radius;
	if (this->NumberOfFilters>0)
		for (i=0;i<this->NumberOfFilters;i++)
			this->Filters[i]->SetLiftingRadius(radius);
}

int vtkMultiresolutionIO::GetLifting()
{
	return this->Lifting;
}

void vtkMultiresolutionIO::InsertNextFilter(vtkWaveletSubdivisionFilter *Filter)
{
	int i;

	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		this->Filters[i+1]=this->Filters[i];
		this->SynthesisMeshes[i+2]=this->Filters[i+1]->GetSubdivisionInput();
	}

	this->SynthesisMeshes[1]=Filter->GetSubdivisionInput();
	this->Filters[0]=Filter;
	this->SynthesisMeshes[0]=Filter->GetOutput();

	this->NumberOfFilters++;
	Filter->SetPointsIds(this->PointsIds);
}

void vtkMultiresolutionIO::EncodeSignificance(int S,int Coord)
{
	this->ArithmeticCoder->EncodeBit(S!=0);
	return;
	this->ArithmeticCoder->EncodeSignificance(S,Coord,0);
	this->SignificanceCodes++;
};

void vtkMultiresolutionIO::EncodeCoeffRefinement(int Coeff)
{
	this->ArithmeticCoder->EncodeSignOrRefinement(Coeff);
};

void vtkMultiresolutionIO::EncodeSign(int Sign)
{
	this->ArithmeticCoder->EncodeSignOrRefinement(Sign);
};

int vtkMultiresolutionIO::DecodeSignificance(int Coord)
{
	return (this->ArithmeticCoder->DecodeBit());
	int S; 
	S=this->ArithmeticCoder->DecodeSignificance(Coord,0);
	this->SignificanceCodes++;
	return (S);
};

int vtkMultiresolutionIO::DecodeCoeffRefinement()
{
	int Coeff; 
	Coeff=this->ArithmeticCoder->DecodeSignOrRefinement();
	return (Coeff);
};

int vtkMultiresolutionIO::DecodeSign()
{
	int Sign;
	Sign=this->ArithmeticCoder->DecodeSignOrRefinement();
	return (Sign);
};


void vtkMultiresolutionIO::GetFirstGoodEdge(vtkIdType Edge,vtkIdType FilterId, vtkIdType &Edge2, vtkIdType &FilterId2, vtkIdType &Vertex)
{

	Edge2=Edge;
	FilterId2=FilterId;

	Vertex=this->EdgeMidPoints[FilterId2]->GetId(Edge);
	if (Vertex>=0)
	{
		return;
	}
	while ((Vertex<0)&&(FilterId2>0))
	{
		Edge2=this->Filters[FilterId2]->TreeFirstChildEdge->GetId(Edge2);
		FilterId2--;
		Vertex=this->EdgeMidPoints[FilterId2]->GetId(Edge2);		
	}
	return;
}

void vtkMultiresolutionIO::ComputeWaveletTree ()
{
	int i,j,k,NumberOfEdges;
	double *Wavelet,*Max,*Max2;
	double AbsoluteWavelet;
	int NumberOfVertices;
	vtkIdType Vertex;
	for (i=0;i<this->NumberOfFilters;i++)
	{
		this->Filters[i]->TreeFirstChildEdge=this->Filters[i]->TreeFirstChildEdge;
		this->TreeNextChildEdge[i]=this->Filters[i]->TreeNextChildEdge;
		this->TreeParentEdge[i]=this->Filters[i]->TreeParentEdge;
		this->Wavelets[i]=this->Filters[i]->Wavelets;
		this->EdgeMidPoints[i]=this->Filters[i]->EdgeMidPoints;

		NumberOfEdges=this->Filters[i]->GetSubdivisionInput()->GetNumberOfEdges();
		this->Maxima[i]=vtkDoubleArray::New();
		this->Maxima[i]->SetNumberOfValues(NumberOfEdges*3);
		for (j=0;j<NumberOfEdges;j++)
		{
			Max=this->Maxima[i]->GetPointer(j*3);
			Max[0]=0;
			Max[1]=0;
			Max[2]=0;
		}
	}

	for (i=0;i<this->NumberOfFilters;i++)
	{
		NumberOfVertices=this->Filters[i]->GetSubdivisionInput()->GetNumberOfPoints();
		NumberOfEdges=this->Filters[i]->GetSubdivisionInput()->GetNumberOfEdges();


		if (i<this->NumberOfFilters-1)
		{
			for (j=0;j<NumberOfEdges;j++)
			{
				Vertex=this->EdgeMidPoints[i]->GetId(j)-NumberOfVertices;
				if (Vertex>=0)
				{
					Max=this->Maxima[i+1]->GetPointer(this->TreeParentEdge[i+1]->GetId(j)*3);
					Wavelet=this->Wavelets[i]->GetPointer(Vertex*3);
					for (k=0;k<3;k++)
					{
						AbsoluteWavelet=Wavelet[k];
						if (fabs(Max[k])<fabs(AbsoluteWavelet))
							Max[k]=AbsoluteWavelet;
					}
				}
			}
		}
		if (i>0)
		{
			NumberOfEdges=this->Filters[i]->GetOutput()->GetNumberOfEdges();	

			for (j=0;j<NumberOfEdges;j++)
			{
				Max=this->Maxima[i]->GetPointer(this->TreeParentEdge[i]->GetId(j)*3);
				Max2=this->Maxima[i-1]->GetPointer(j*3);
				for (k=0;k<3;k++)
				{
					if (fabs(Max[k])<fabs(Max2[k]))
						Max[k]=Max2[k];
				}
			}
		}
	}
}

void vtkMultiresolutionIO::EncodeProgressivePrecision()
{
	this->ComputeWaveletTree();

	vtkIdType i,j,BitPlane;
	WaveletCoefficientSet Set;
	WaveletCoefficient Coeff;

	double Max[3],*Wav,AbsoluteCoeff;
	vtkIdType Vertex,Edge;
	vtkIdType NumberOfVertices,FilterId;
	int Count=0; 

	for (i=0;i<3;i++)
	{
		this->InsignificantCoeffs[i].clear();
		this->InsignificantSets[i].clear();
		this->SignificantCoeffs[i].clear();
	}
	this->ArithmeticCoder->StartCoding();
	this->ArithmeticCoder->EncodeByte((this->NumberOfBitPlanes<<4)+this->NumberOfStartBitPlanes);
	this->ArithmeticCoder->InitConnectivityQSModels(1);

	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		if (this->Filters[i]->GetSubdivisionType()==0)
		{
			this->ArithmeticCoder->EncodeBit(0);
		}
		else
		{
			if (this->Filters[i]->GetSubdivisionType()==1)
			{
				this->ArithmeticCoder->EncodeBit(1);
				this->ArithmeticCoder->EncodeBit(0);
			}
			else
			{
				this->ArithmeticCoder->EncodeBit(1);
				this->ArithmeticCoder->EncodeBit(1);
			}
		}

		this->Filters[i]->SetIOType(1);
		this->Filters[i]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[i]->Subdivide();
	}
	
	this->ConnectivityBytes=this->ArithmeticCoder->StopCoding();
	this->ArithmeticCoder->StartCoding();

	Max[0]=0;
	Max[1]=0;
	Max[2]=0;

	NumberOfVertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();

	for (i=0;i<this->SynthesisMeshes[this->NumberOfFilters]->GetNumberOfEdges();i++)
	{
		Wav=this->Maxima[this->NumberOfFilters-1]->GetPointer(i*3);
		for (j=0;j<3;j++)
		{
			AbsoluteCoeff=Wav[j];
			if (fabs(Max[j])<fabs(AbsoluteCoeff))
				Max[j]=AbsoluteCoeff;
		}
		this->GetFirstGoodEdge(i,this->NumberOfFilters-1,Coeff.Edge,Coeff.FilterId,Vertex);

		if (Vertex>=0)
		{
			Vertex-=this->Filters[Coeff.FilterId]->GetSubdivisionInput()->GetNumberOfPoints();
			Coeff.WaveletId=Vertex;
			// Add entry to LIP
			Wav=this->Wavelets[Coeff.FilterId]->GetPointer(Vertex*3);
			for (j=0;j<3;j++)
			{
				AbsoluteCoeff=Wav[j];
				if (fabs(Max[j])<fabs(AbsoluteCoeff))
					Max[j]=AbsoluteCoeff;
				Coeff.Wavelet=AbsoluteCoeff;
				Coeff.Count=Count++;
				Coeff.Sn=2;
				InsignificantCoeffs[j].push_back(Coeff);
			}

			// Add entry to LIS
			Set.Edge=Coeff.Edge;
			Set.FilterId=Coeff.FilterId;
			Set.Type=1;
			Set.Sn=2;
			Set.Count=Count++;

			Wav=this->Maxima[Set.FilterId]->GetPointer(Set.Edge*3);
			if (Set.FilterId>0)
				for (j=0;j<3;j++)
				{
					Set.Max=Wav[j];
					InsignificantSets[j].push_back(Set);
				}
		}
	}

	AbsoluteCoeff=Max[0];
	if (fabs(AbsoluteCoeff)<fabs(Max[1]))
		AbsoluteCoeff=Max[1];
	if (fabs(AbsoluteCoeff)<fabs(Max[2]))
		AbsoluteCoeff=Max[2];

	int N,Descendants,Sn;
	vtkIdType Edge2,FilterId2;
	double T,MaxCoeff,TestCoeff;
	std::list<WaveletCoefficientSet>::iterator SetIter;
	std::list<WaveletCoefficient>::iterator CoeffIter;
	std::list<WaveletCoefficientSet>::iterator SetIter2;
	std::list<WaveletCoefficient>::iterator CoeffIter2;

	N=(int) floor(log(fabs(AbsoluteCoeff))/log(2.00));
	for (i=0;i<3;i++)
	{
		for (CoeffIter=InsignificantCoeffs[i].begin();
			CoeffIter!=InsignificantCoeffs[i].end();CoeffIter++)
			CoeffIter->N=N;
	}

	for (i=0;i<this->NumberOfFilters;i++)
	{
		this->SauvWavelets[i]=vtkDoubleArray::New();
		this->SauvWavelets[i]->DeepCopy(Wavelets[i]);
	}

	this->ArithmeticCoder->EncodeByte(N+128);
	this->ArithmeticCoder->EncodeByte(this->NumberOfBitPlanes);
	this->ArithmeticCoder->EncodeByte(this->NumberOfStartBitPlanes);
	T=pow(2.0,N);

	for (BitPlane=0;BitPlane<this->NumberOfBitPlanes;BitPlane++)
	{
		if (BitPlane>=this->NumberOfStartBitPlanes)
			this->ArithmeticCoder->StartCoding();
		this->ArithmeticCoder->InitZerotreeQSModels(1);
		for (i=0;i<3;i++)
		{
			int number=0;

			//2.1 For Each Entry (i,j) in the LIP
			for (CoeffIter=InsignificantCoeffs[i].begin();
				CoeffIter!=InsignificantCoeffs[i].end();)
			{
				number++;

				Sn=(fabs(CoeffIter->Wavelet)>=T);
				this->EncodeSignificance(Sn,i);

				//2.1.1 Output Sn(i,j)
				if (Sn==0)
					CoeffIter++;
				else
				{
					// If Sn(i,j)=1 then move (i,j) to the LSP and output the sign of cij
					Sn=(CoeffIter->Wavelet<0);
					this->EncodeSign(Sn);

					Coeff.Edge=CoeffIter->Edge;
					Coeff.FilterId=CoeffIter->FilterId;
					Coeff.Wavelet=CoeffIter->Wavelet;
					Coeff.WaveletId=CoeffIter->WaveletId;
					Coeff.N=N;
					SignificantCoeffs[i].push_back(Coeff);

					CoeffIter=InsignificantCoeffs[i].erase(CoeffIter);
				}
			}

			//2.2 for each entry (i,j) in the LIS do:
			for (SetIter=InsignificantSets[i].begin();
				SetIter!=InsignificantSets[i].end();)
			{
				if (SetIter->Type==1)
				{
					// 2.2.1 if the entry is of type A
					// Output Sn(D(i,j))

					Sn=(fabs(SetIter->Max)>=T);
					this->EncodeSignificance(Sn,i);

					if (Sn==0)
						SetIter++;
					else
					{
						// is Sn(i,j)=1 then:
						Descendants=0;
						MaxCoeff=0;
						FilterId=SetIter->FilterId-1;

						// for each (k,l) in O(i,j)
						Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
						while (Edge>=0)
						{
							this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
							TestCoeff=this->Maxima[FilterId2]->GetPointer(Edge2*3)[i];
							if (fabs(MaxCoeff)<fabs(TestCoeff))
								MaxCoeff=TestCoeff;

							if (Vertex>=0)
							{
								if (FilterId2>0)
									Descendants++;
								Coeff.WaveletId=Vertex-this->Filters[FilterId2]->GetSubdivisionInput()->GetNumberOfPoints();
								Coeff.Edge=Edge2;
								Coeff.FilterId=FilterId2;
								Coeff.Wavelet=this->Wavelets[FilterId2]->GetPointer(Coeff.WaveletId*3)[i];
								Coeff.N=N;
								Coeff.Count=Count;
								Coeff.Sn=2;

								//Output Sn(k,l)
								Sn=(fabs(Coeff.Wavelet)>=T);

								this->EncodeSignificance(Sn,i);
								if (Sn==0)
									// if Sn(k,l)=0 then add (k,l) to the end of the LIP
									InsignificantCoeffs[i].push_back(Coeff);
								else
								{
									// if Sn(k,l)=1 then add (k,l) to the end of the LSP and output the sign of cij
									Sn=Coeff.Wavelet<0;
									this->EncodeSign(Sn);
									SignificantCoeffs[i].push_back(Coeff);
								}

							}
							if (FilterId>=0)
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							else
								Edge=-1;
						}
						Count++;
						// if L(i,j) non empty
						if ((Descendants>0)&&(FilterId>0))
						{

							SetIter->Max=MaxCoeff;
							SetIter->Type=0;
						}
						else
						{
							// else remove (i,j) from the LIS
							SetIter=InsignificantSets[i].erase(SetIter);
						}
					}
				}

				// 2.2.2 if the entry is of type B then
				if (SetIter!=InsignificantSets[i].end())
				{
					if (SetIter->Type==0)
					{
						Sn=(fabs(SetIter->Max)>=T);
						// Output Sn(L(i,j))
						this->EncodeSignificance(Sn,i);

						if (Sn==0)
							SetIter++;
						else
						{
							// if Sn(L(i,j))=1 then
							MaxCoeff=0;
							Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
							FilterId=SetIter->FilterId-1;
							// add each (k,l) In O(i,j) to the end of the LIS as an entry of type A
							while (Edge>=0)
							{
								this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
								if ((Vertex>=0)&&(FilterId2>0))
								{
									Set.Edge=Edge2;
									Set.FilterId=FilterId2;
									Set.Max=this->Maxima[FilterId2]->GetPointer(Edge2*3)[i];
									Set.Type=1;
									Set.Sn=2;
									Set.Count=Count;
									InsignificantSets[i].push_back(Set);
								}
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							}
							// remove (i,j) from the LIS;
							SetIter=InsignificantSets[i].erase(SetIter);
							Count++;
						}
					}
				}
			}		

			// Refinement pass
			int newnumber=0,k,Ptest;
			number=0;
			for(CoeffIter=SignificantCoeffs[i].begin();
				CoeffIter!=SignificantCoeffs[i].end();CoeffIter++)
			{
				number++;
				if (CoeffIter->N!=N)
				{
					newnumber++;
					TestCoeff=floor(fabs(CoeffIter->Wavelet)/T);
					Ptest=(int) TestCoeff;
					if (Ptest%2==1)
					{
						this->EncodeCoeffRefinement(1);
					}
					else
					{
						this->EncodeCoeffRefinement(0);
					}

				}
			}

			number=0;
			for (k=0;k<this->NumberOfFilters;k++)
			{
				for (j=0;j<this->Filters[k]->GetOutput()->GetNumberOfPoints()
					-this->Filters[k]->GetSubdivisionInput()->GetNumberOfPoints();j++)
				{
					Wav=this->Wavelets[k]->GetPointer(j*3);
					if (fabs(Wav[i])>=T)
						number++;
				}
			}

		}
		N--;
		T=T/2;
		if (BitPlane>=this->NumberOfStartBitPlanes-1)
			this->DataSize[BitPlane]=this->ArithmeticCoder->StopCoding();
		else
			this->DataSize[BitPlane]=0;
	}

	this->ArithmeticCoder->CloseFile();
	if (this->WriteRepport > 0)
	{
		std::string outinfo = this->InputFileName;
		if (this->Lifting > 0)
			outinfo = "lifting_repport.txt";
		else
			outinfo = "repport.txt";

		std::ofstream Repport;
		Repport.open(outinfo, ofstream::out | ofstream::trunc);

		double duration,ratio;
#if ( (VTK_MAJOR_VERSION >= 5))
		duration = this->Timer->GetUniversalTime()- this->StartTime;
#else
		duration = this->Timer->GetCurrentTime()- this->StartTime;
#endif
		ratio= ((double)this->Filters[0]->GetOutput()->GetNumberOfCells())/duration;

		int Sum;
		int NumberOfVertices=this->SynthesisMeshes[0]->GetNumberOfPoints();
		Repport<<"Zerotree Encoding of a mesh with "<<this->SynthesisMeshes[0]->GetNumberOfCells()<<
			" faces and "<<NumberOfVertices<<" vertices, "<<this->NumberOfFilters<<" resolution levels"<<endl;

		Repport<<"Coding of mesh connectivity+ base mesh geometry : "<<this->ConnectivityBytes<<
			" bytes ( "<<this->ConnectivityBytes*8<<" bits, "<<
			8.0*(double)this->ConnectivityBytes/(double)NumberOfVertices<<" bits/vertex)"<<endl;
		Repport<<"Encoding time : "<<duration<<" seconds ( "<<ratio<<" faces/s)"<<endl;

		Sum=this->ConnectivityBytes+this->DataSize[this->NumberOfStartBitPlanes-1];

		Repport<<"Bitplanes 0 to "<<this->NumberOfStartBitPlanes-1<<" : "
			<<this->DataSize[this->NumberOfStartBitPlanes-1]<<" Bytes, Total Data : "<<
			Sum<<" Bytes ( "<<Sum*8<<" bits, "<<8.0*(double)Sum/(double)NumberOfVertices<<" bits/Vertex)"<<endl;

		for (BitPlane=this->NumberOfStartBitPlanes;BitPlane<this->NumberOfBitPlanes;BitPlane++)
		{
			Sum+=this->DataSize[BitPlane];
			Repport<<"Bitplane "<<BitPlane<<" : "<<this->DataSize[BitPlane]<<" Bytes, Total Data : "<<
				Sum<<" Bytes ( "<<Sum*8<<" bits, "<<8.0*(double)Sum/(double)NumberOfVertices<<" bits/Vertex)"<<endl;
		}
		Repport.close();
	}
}

void vtkMultiresolutionIO::DecodeProgressivePrecision()
{
	int i,j,BitPlane;
	WaveletCoefficientSet Set;
	WaveletCoefficient Coeff;

	vtkIdType Vertex,Edge;
	int NumberOfVertices,FilterId,Count=0;

	for (i=0;i<3;i++)
	{
		this->InsignificantCoeffs[i].clear();
		this->InsignificantSets[i].clear();
		this->SignificantCoeffs[i].clear();
	}

	this->ArithmeticCoder->StartDecoding();

	int SubdivisionType;
	int BitPlanesCode=this->ArithmeticCoder->DecodeByte();
	this->NumberOfBitPlanes=BitPlanesCode>>4;
	this->NumberOfStartBitPlanes=BitPlanesCode-(this->NumberOfBitPlanes<<4);

	this->ArithmeticCoder->InitConnectivityQSModels(0);

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		if (this->ArithmeticCoder->DecodeBit()==0)
			SubdivisionType=0;
		else
		{
			if (this->ArithmeticCoder->DecodeBit()==0)
				SubdivisionType=1;
			else
				SubdivisionType=2;
		}
		
		this->Filters[j]=this->NewFilter(SubdivisionType);

		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);

		this->Filters[j]->SetSubdivisionType(SubdivisionType);
		this->Filters[j]->SetIOType(2);	
		this->Filters[j]->SetLifting(this->Lifting);
		this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[j]->GeometryPrediction=this->GeometryPrediction;
		this->Filters[j]->SetArithmeticType(0);
		this->EdgeMidPoints[j]=this->Filters[j]->EdgeMidPoints;
		this->TreeNextChildEdge[j]=this->Filters[j]->TreeNextChildEdge;

		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->Subdivide();
		this->Filters[j]->Wavelets=vtkDoubleArray::New();
		this->Filters[j]->Wavelets->SetNumberOfComponents(3);
		this->Filters[j]->Wavelets->SetNumberOfTuples(this->Filters[j]->GetOutput()->GetNumberOfPoints()-
			this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints());

		this->Wavelets[j]=this->Filters[j]->Wavelets;
		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();

		if (this->DisplayText)
			cout<<"Level "<<this->NumberOfFilters-j<<": "<<this->SynthesisMeshes[j]->GetNumberOfCells()
			<<" faces, "<<this->SynthesisMeshes[j]->GetNumberOfPoints()<<" vertices "<<endl;

	}

	int k;
	double Wav[3];
	for (k=0;k<this->NumberOfFilters;k++)
	{
		for (j=0;j<this->Filters[k]->GetOutput()->GetNumberOfPoints()
			-this->Filters[k]->GetSubdivisionInput()->GetNumberOfPoints();j++)
		{
			this->Wavelets[k]->GetTuple(j,Wav);
			Wav[0]=0;
			Wav[1]=0;
			Wav[2]=0;
			this->Wavelets[k]->SetTuple(j,Wav);
		}
	}
	this->Reconstruct();
	if (this->DisplayText)
		cout<<"Final mesh: "<<this->SynthesisMeshes[0]->GetNumberOfCells()<<" faces"<<endl;
	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
		this->MeshWindow->Render();
		this->MeshWindow->SetWindowName("Progressive precision reconstruction");
		cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
		this->MeshWindow->Interact();

	}
	this->ArithmeticCoder->StopDecoding();

	this->ArithmeticCoder->StartDecoding();

	NumberOfVertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();
	for (i=0;i<this->SynthesisMeshes[this->NumberOfFilters]->GetNumberOfEdges();i++)
	{
		this->GetFirstGoodEdge(i,this->NumberOfFilters-1,Coeff.Edge,Coeff.FilterId,Vertex);
		if (Vertex>=0)
		{
			Vertex-=this->Filters[Coeff.FilterId]->GetSubdivisionInput()->GetNumberOfPoints();
			Coeff.WaveletId=Vertex;
			Coeff.Count=Count++;
			Coeff.Sn=2;
			// Add entry to LIP
			for (j=0;j<3;j++)
				InsignificantCoeffs[j].push_back(Coeff);

			// Add entry to LIS
			Set.Edge=Coeff.Edge;
			Set.FilterId=Coeff.FilterId;
			Set.Type=1;
			Set.Count=Count++;
			Set.Sn=2;

			if (Set.FilterId>0)
				for (j=0;j<3;j++)
					InsignificantSets[j].push_back(Set);
		}
		else
		{
//			cout<<"Bizarre!!!"<<endl;
		}

	}


	int N,Descendants,Sn;
	vtkIdType Edge2,FilterId2;
	double T,MaxCoeff;
	std::list<WaveletCoefficientSet>::iterator SetIter;
	std::list<WaveletCoefficient>::iterator CoeffIter;

	N=this->ArithmeticCoder->DecodeByte()-128;
	this->NumberOfBitPlanes=this->ArithmeticCoder->DecodeByte();
	this->NumberOfStartBitPlanes=this->ArithmeticCoder->DecodeByte();

	T=pow(2.0,N);
	for (i=0;i<3;i++)
	{
		for (CoeffIter=InsignificantCoeffs[i].begin();
			CoeffIter!=InsignificantCoeffs[i].end();CoeffIter++)
			CoeffIter->N=N;
	}


	for (BitPlane=0;BitPlane<this->NumberOfBitPlanes;BitPlane++)
	{
		if (BitPlane>=this->NumberOfStartBitPlanes)
			this->ArithmeticCoder->StartDecoding();
		this->ArithmeticCoder->InitZerotreeQSModels(0);

		for (i=0;i<3;i++)
		{
			int number=0;
			//2.1 For Each Entry (i,j) in the LIP
			for (CoeffIter=InsignificantCoeffs[i].begin();
				CoeffIter!=InsignificantCoeffs[i].end();)
			{
				number++;
				//2.1.1 Input Sn(i,j)
				Sn=this->DecodeSignificance(i);

				if (Sn==0)
					CoeffIter++;
				else
				{

					// If Sn(i,j)=1 then move (i,j) to the LSP and input the sign of cij
					Sn=this->DecodeSign();
					if (Sn==1)
						this->Wavelets[CoeffIter->FilterId]->GetPointer(CoeffIter->WaveletId*3)[i]=-1.5*pow(2.0,N);
					else
						this->Wavelets[CoeffIter->FilterId]->GetPointer(CoeffIter->WaveletId*3)[i]=1.5*pow(2.0,N);

					Coeff.WaveletId=CoeffIter->WaveletId;
					Coeff.Edge=CoeffIter->Edge;
					Coeff.FilterId=CoeffIter->FilterId;
					Coeff.N=N;
					SignificantCoeffs[i].push_back(Coeff);
					CoeffIter=InsignificantCoeffs[i].erase(CoeffIter);
				}
			}

			//2.2 for each entry (i,j) in the LIS do:
			for (SetIter=InsignificantSets[i].begin();
				SetIter!=InsignificantSets[i].end();)
			{
				if (SetIter->Type==1)
				{
					// 2.2.1 if the entry is of type A
					// Output Sn(D(i,j))
					Sn=this->DecodeSignificance(i);


					if (Sn==0)
						SetIter++;
					else
					{
						// is Sn(i,j)=1 then:
						Descendants=0;
						MaxCoeff=0;
						FilterId=SetIter->FilterId-1;


						// for each (k,l) in O(i,j)
						Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
						while (Edge>=0)
						{
							this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
							if (Vertex>=0)
							{
								if (FilterId2>0)
									Descendants++;

								Coeff.WaveletId=Vertex-this->Filters[FilterId2]->GetSubdivisionInput()->GetNumberOfPoints();
								Coeff.Edge=Edge2;
								Coeff.FilterId=FilterId2;
								Coeff.N=N;


								//Input Sn(k,l)
								Sn=this->DecodeSignificance(i);

								if (Sn==0)
								{
									// if Sn(k,l)=0 then add (k,l) to the end of the LIP
									InsignificantCoeffs[i].push_back(Coeff);
								}
								else
								{
									// if Sn(k,l)=1 then add (k,l) to the end of the LSP and output the sign of cij
									Sn=this->DecodeSign();
									if (Sn==1)
										this->Wavelets[Coeff.FilterId]->GetPointer(Coeff.WaveletId*3)[i]=-1.5*pow(2.0,N);
									else
										this->Wavelets[Coeff.FilterId]->GetPointer(Coeff.WaveletId*3)[i]=1.5*pow(2.0,N);

									SignificantCoeffs[i].push_back(Coeff);
								}
							}
							if (FilterId>=0)
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							else
								Edge=-1;
						}
						// if L(i,j) non empty
						if ((Descendants>0)&&(FilterId>0))
						{
							SetIter->Type=0;
						}
						else
						{
							// else remove (i,j) from the LIS
							SetIter=InsignificantSets[i].erase(SetIter);
						}
					}
				}

				// 2.2.2 if the entry is of type B then
				if (SetIter!=InsignificantSets[i].end())
				{
					if (SetIter->Type==0)
					{
						Sn=this->DecodeSignificance(i);

						// Input Sn(L(i,j))
						if (Sn==0)
							SetIter++;
						else
						{
							// if Sn(L(i,j))=1 then
							MaxCoeff=0;
							Edge=this->Filters[SetIter->FilterId]->TreeFirstChildEdge->GetId(SetIter->Edge);
							FilterId=SetIter->FilterId-1;
							// add each (k,l) In O(i,j) to the end of the LIS as an entry of type A
							while (Edge>=0)
							{
								this->GetFirstGoodEdge(Edge,FilterId,Edge2,FilterId2,Vertex);
								if ((Vertex>=0)&&(FilterId2>0))
								{
									Set.Edge=Edge2;
									Set.FilterId=FilterId2;
									Set.Type=1;
									InsignificantSets[i].push_back(Set);
								}
								Edge=this->TreeNextChildEdge[SetIter->FilterId]->GetId(Edge);
							}
							// remove (i,j) from the LIS;
							SetIter=InsignificantSets[i].erase(SetIter);
						}
					}
				}
			}		

			// Refinement pass
			int newnumber=0,ref;
			double *WaveletPointer;
			number=0;
			for(CoeffIter=SignificantCoeffs[i].begin();
				CoeffIter!=SignificantCoeffs[i].end();CoeffIter++)
			{
				number++;
				if (CoeffIter->N!=N)
				{
					newnumber++;
					WaveletPointer=this->Wavelets[CoeffIter->FilterId]->GetPointer(CoeffIter->WaveletId*3);
					ref=this->DecodeCoeffRefinement();

					if (WaveletPointer[i]<0)
					{
						if (ref==1)
							WaveletPointer[i]-=pow(2.0, N-1);
						else
							WaveletPointer[i]+=pow(2.0, N-1);
					}
					else
					{
						if (ref==1)
							WaveletPointer[i]+=pow(2.0, N-1);
						else
							WaveletPointer[i]-=pow(2.0,N-1);
					}
				}
			}
		}

		N--;

		cout<<"Bitplane :"<<BitPlane<<", Threshold="<<T<<endl;
		T=T/2;

		this->SynthesisMeshes[0]->Modified();
		if (BitPlane>=this->NumberOfStartBitPlanes-1)
		{
			this->ArithmeticCoder->StopDecoding();
			this->Reconstruct();
			if (this->Display!=0)
			{
				this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
				this->MeshWindow->Render();
			}
			if ((this->Display==2)||((BitPlane==this->NumberOfStartBitPlanes-1)&&(this->Display!=0)))
				this->MeshWindow->Interact();
			if (this->WriteOutput==1)
			{
				if (this->FileType==0)
				{
					std::stringstream strfile;
					strfile<<"Mesh"<<BitPlane<<".iv";
					this->Filters[0]->GetOutput()->WriteInventor(strfile.str().c_str());
				}
				else
				{
					std::stringstream strfile;
					strfile<<"Mesh"<<BitPlane<<".ply";
					vtkPLYWriter *Writer=vtkPLYWriter::New();
					Writer->SetInputData(this->Filters[0]->GetOutput());
					Writer->SetFileName(strfile.str().c_str());
					//Writer->SetFileTypeToASCII();
					Writer->Write();
					Writer->Delete();
				}
			}
			if (this->Capture==1)
			{
				std::stringstream strfile;
				strfile<<"Mesh"<<BitPlane<<".bmp";
				this->MeshWindow->Capture(strfile.str().c_str());
			}
		}
	}

	if (this->Display!=0)
	{
		this->MeshWindow->Render();
		this->MeshWindow->Interact();
	}
	this->ArithmeticCoder->CloseFile();
	this->Output=this->SynthesisMeshes[0];
}

void vtkMultiresolutionIO::Execute()
{

}

void vtkMultiresolutionIO::Approximate()
{
	int i;
	for (i=0;i<this->NumberOfFilters;i++)
	{
		this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
		this->Filters[i]->SetLifting(this->Lifting);
		this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[i]->Approximate();

	}
}

void vtkMultiresolutionIO::Reconstruct()
{
	int i;
	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
		this->Filters[i]->SetLifting(this->Lifting);
		this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[i]->Reconstruct();
		this->SynthesisMeshes[i]=this->Filters[i]->GetOutput();
	}
	this->SynthesisMeshes[0]->Modified();
	this->SynthesisMeshes[0]->GetPoints()->Modified();
	this->Output=this->SynthesisMeshes[0];
	this->Output->Register(this);
}
void vtkMultiresolutionIO::Read()
{
	this->ArithmeticCoder=vtkArithmeticCoder::New();
	this->ArithmeticCoder->OpenFile(this->FileName,0);
	this->ArithmeticCoder->InitConnectivityQSModels(0);
	this->ArithmeticCoder->StartDecoding();
	this->ArithmeticType=this->ArithmeticCoder->DecodeBit();

	vtkSurface *BaseMesh=this->DecodeMeshConnectivity();
	this->DecodeMeshGeometry(BaseMesh);
	this->DecodeLiftingProperties();
	this->NumberOfFilters=this->ArithmeticCoder->DecodeByte();

	int Type,i;
	int bit1,bit2;
	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		bit1=this->ArithmeticCoder->DecodeBit();
		bit2=this->ArithmeticCoder->DecodeBit();
		Type=bit1*2+bit2;

		this->Filters[i]=this->NewFilter(Type);
		this->Filters[i]->SetSubdivisionType(Type);
	}


	this->ArithmeticCoder->StopDecoding();


	if (this->DisplayText)
		cout<<"Level 0 : "<<BaseMesh->GetNumberOfCells()
		<<" faces, "<<BaseMesh->GetNumberOfPoints()<<" vertices and "
		<<BaseMesh->GetNumberOfEdges()<<" edges "<<endl;

	this->SynthesisMeshes[this->NumberOfFilters]=BaseMesh;
	if (this->ArithmeticType!=0)
        this->DecodeProgressiveResolution();
	else
		this->DecodeProgressivePrecision();
}

void vtkMultiresolutionIO::ReadPerMesh()
{
	std::string outinfo = this->InputFileName;
	clock_t start, end;

	this->ArithmeticCoder = vtkArithmeticCoder::New();
	this->ArithmeticCoder->OpenFile(this->FileName, 0);
	this->ArithmeticCoder->InitConnectivityQSModels(0);

	start = clock();

	this->ArithmeticCoder->StartDecoding();
	this->ArithmeticType = this->ArithmeticCoder->DecodeBit();

	vtkSurface* BaseMesh = this->DecodeMeshConnectivity();
	this->DecodeMeshGeometry(BaseMesh);
	this->DecodeLiftingProperties();
	this->NumberOfFilters = this->ArithmeticCoder->DecodeByte();

	int Type, i;
	int bit1, bit2;
	for (i = this->NumberOfFilters - 1; i >= 0; i--)
	{
		bit1 = this->ArithmeticCoder->DecodeBit();
		bit2 = this->ArithmeticCoder->DecodeBit();
		Type = bit1 * 2 + bit2;

		this->Filters[i] = this->NewFilter(Type);
		this->Filters[i]->SetSubdivisionType(Type);
	}

	this->ArithmeticCoder->StopDecoding();

	end = clock();
	double reading_time = (double)(end - start) / CLOCKS_PER_SEC;

	start = clock();
	this->GetRenderWindow()->Render();
	end = clock();
	double rendering_time = (double)(end - start) / CLOCKS_PER_SEC;

	// file name setting
	if (this->Lifting > 0)
		outinfo += "_DS_lift_0.txt";
	else
		outinfo += "_DS_nolift.txt";

	std::ofstream out;

	// for each resolution in multi-resolution mesh:
	//   (1) print reading time, rendering time, (reading + rendering) time
	out.open(outinfo, ios::app);
	out << "level \t" << "reading_time  \t" << "rendering_time \t" << "total_time" << endl;
	out << this->NumberOfFilters << " \t" << reading_time << " \t" << rendering_time << " \t" << reading_time + rendering_time << endl;
	out.close();

	this->SynthesisMeshes[this->NumberOfFilters] = BaseMesh;
	if (this->ArithmeticType != 0)
		this->DecodeProgressiveResolutionDS();
	else
		this->DecodeProgressivePrecision();
}

void vtkMultiresolutionIO::DecodeProgressiveResolution()
{
	int j;
	double start2,finish2;

	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[this->NumberOfFilters]);
		this->MeshWindow->Render();
		this->MeshWindow->SetWindowName("Progressive resolution reconstruction");
#if ( (VTK_MAJOR_VERSION >= 5))
		start2=this->Timer->GetUniversalTime();
		finish2=this->Timer->GetUniversalTime();
		while (finish2 - start2<this->Time)
			finish2=this->Timer->GetUniversalTime();

#else
		start2=this->Timer->GetCurrentTime();
		finish2=this->Timer->GetCurrentTime();
		while (finish2 - start2<this->Time)
			finish2=this->Timer->GetCurrentTime();
#endif			
		this->MeshWindow->Interact();
	}

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);
		this->Filters[j]->SetIOType(2);	
		this->Filters[j]->SetLifting(this->Lifting);
		this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[j]->SetArithmeticType(1);
		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->Quantization=this->Quantization;

		this->ArithmeticCoder->StartDecoding();
		if (j==this->NumberOfFilters-1)
			this->DecodeScalingFactors(this->Filters[j]->GetSubdivisionInput());


		this->Filters[j]->Subdivide();
		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();
		this->ArithmeticCoder->StopDecoding();
		this->ArithmeticCoder->StartDecoding();
		this->Filters[j]->ReadCoefficients();
		this->ArithmeticCoder->StopDecoding();
		this->Filters[j]->Reconstruct();

		if (this->DisplayText)
			cout<<"Level "<<this->NumberOfFilters-j<<": "<<this->SynthesisMeshes[j]->GetNumberOfCells()
			<<" faces, "<<this->SynthesisMeshes[j]->GetNumberOfPoints()<<" vertices and "
			<<this->SynthesisMeshes[j]->GetNumberOfEdges()<<" edges "<<endl;


		if (this->Display!=0)
		{
			this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
			this->MeshWindow->Render();

#if ( (VTK_MAJOR_VERSION >= 5))
			start2=this->Timer->GetUniversalTime();
			finish2=this->Timer->GetUniversalTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetUniversalTime();

#else
			start2=this->Timer->GetCurrentTime();
			finish2=this->Timer->GetCurrentTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetCurrentTime();
#endif
			if (this->Capture==1)
			{
				std::stringstream strfile;
				strfile<<"Mesh"<<this->NumberOfFilters-j<<".bmp";
				this->MeshWindow->Capture(strfile.str().c_str());
			}

			if (this->Display==2)
			{
				cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
				this->MeshWindow->Interact();
			}
		}

	}

	if (this->DisplayText)
		cout<<"Final mesh: "<<this->SynthesisMeshes[0]->GetNumberOfCells()<<" faces"<<endl;
	this->ArithmeticCoder->CloseFile();

	this->SynthesisMeshes[0]->Modified();
	if (this->Display!=0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
		this->MeshWindow->Render();
		cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
		this->MeshWindow->Interact();
	}
	if (this->WriteOutput==-1)
	{
		// output file name initialization
		std::string strout;
		if (this->InputFileName != NULL)
		{
			strout = this->InputFileName;
			strout += "_";
		}
		else
			strout = "Mesh";

		this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->UnQuantizeCoordinates();
		if (this->FileType == 0)
		{
			std::stringstream strfile;
			strfile << strout << "0.iv";
			this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->WriteInventor(strfile.str().c_str());
		}
		else
		{
			vtkPLYWriter *Writer=vtkPLYWriter::New();
			std::stringstream strfile;
			strfile << strout << "0.ply";
			Writer->SetFileName(strfile.str().c_str());
			Writer->SetInputData(this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput());
			//Writer->SetFileTypeToASCII();
			Writer->Write();
			Writer->Delete();
		}
		double Factor,Tx,Ty,Tz;
		this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetScalingFactors(Factor,Tx,Ty,Tz);
		for (j=0;j<this->NumberOfFilters;j++)
		{
			this->Filters[j]->GetOutput()->SetScalingFactors(Factor,Tx,Ty,Tz);
			this->Filters[j]->GetOutput()->UnQuantizeCoordinates();
			std::stringstream strfile;
			if (this->FileType==0)
			{
				strfile << strout << j+1 << ".iv";
				this->Filters[j]->GetOutput()->WriteInventor(strfile.str().c_str());
			}
			else
			{
				strfile << strout << j+1 << ".ply";
				vtkPLYWriter *Writer=vtkPLYWriter::New();
				Writer->SetFileName(strfile.str().c_str());
				Writer->SetInputData(this->Filters[j]->GetOutput());
				//Writer->SetFileTypeToASCII();
				Writer->Write();
				Writer->Delete();
			}
		}
	}
	if (this->WriteOutput >= 1)
	{
		// output file name initialization
		std::string strout;
		if (this->InputFileName != NULL)
		{
			strout = this->InputFileName;
			strout += "_";
		}
		else
			strout = "Mesh";

		double Factor, Tx, Ty, Tz;
		this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetScalingFactors(Factor, Tx, Ty, Tz);
		if (this->FileType == 0)
		{
			this->Filters[0]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
			this->Filters[0]->GetOutput()->UnQuantizeCoordinates();

			std::stringstream strfile;
			strfile << strout << "0.iv";
			this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->WriteInventor(strfile.str().c_str());
		}
		else
		{
			for (j = 0; j < MIN(this->WriteOutput,this->NumberOfFilters); j++)
			{
				this->Filters[j]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
				this->Filters[j]->GetOutput()->UnQuantizeCoordinates();

				std::stringstream strfile;
				if (this->FileType == 0)
				{
					strfile << strout << j+1 << ".iv";
					this->Filters[j]->GetOutput()->WriteInventor(strfile.str().c_str());
				}
				else
				{
					strfile << strout << j+1 << ".ply";
					vtkPLYWriter* Writer = vtkPLYWriter::New();
					Writer->SetFileName(strfile.str().c_str());
					Writer->SetInputData(this->Filters[j]->GetOutput());
					//Writer->SetFileTypeToASCII();
					Writer->Write();
					Writer->Delete();
				}
			}
		}
	}
}

///////////////////// for DS only //////////////////////
void vtkMultiresolutionIO::DecodeProgressiveResolutionDS()
{
	int j;
	double start2, finish2;

	if (this->Display != 0)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[this->NumberOfFilters]);
		this->MeshWindow->Render();
		this->MeshWindow->SetWindowName("Progressive resolution reconstruction");
#if ( (VTK_MAJOR_VERSION >= 5))
		start2 = this->Timer->GetUniversalTime();
		finish2 = this->Timer->GetUniversalTime();
		while (finish2 - start2 < this->Time)
			finish2 = this->Timer->GetUniversalTime();

#else
		start2 = this->Timer->GetCurrentTime();
		finish2 = this->Timer->GetCurrentTime();
		while (finish2 - start2 < this->Time)
			finish2 = this->Timer->GetCurrentTime();
#endif			
		this->MeshWindow->Interact();
	}

	std::string outinfo = this->InputFileName;
	if (this->Lifting > 0)
		outinfo += "_DS_lift_0.txt";
	else
		outinfo += "_DS_nolift.txt";

	std::ofstream out;

	clock_t start, end;

	for (j = this->NumberOfFilters - 1; j >= 0; j--)
	{
		this->Filters[j]->SetInput(this->SynthesisMeshes[j + 1]);
		this->Filters[j]->SetIOType(2);
		this->Filters[j]->SetLifting(this->Lifting);
		this->Filters[j]->SetLiftingRadius(this->LiftingRadius);
		this->Filters[j]->SetArithmeticType(1);
		this->Filters[j]->ArithmeticCoder = this->ArithmeticCoder;
		this->Filters[j]->Quantization = this->Quantization;

		start = clock();

		this->ArithmeticCoder->StartDecoding();
		if (j == this->NumberOfFilters - 1)
			this->DecodeScalingFactors(this->Filters[j]->GetSubdivisionInput());


		this->Filters[j]->Subdivide();
		this->SynthesisMeshes[j] = this->Filters[j]->GetOutput();
		this->ArithmeticCoder->StopDecoding();
		this->ArithmeticCoder->StartDecoding();
		this->Filters[j]->ReadCoefficients();
		this->ArithmeticCoder->StopDecoding();
		this->Filters[j]->Reconstruct();

		end = clock();
		double reading_time = (double)(end - start) / CLOCKS_PER_SEC;

		start = clock();
		this->MeshWindow->Render();
		end = clock();
		double rendering_time = (double)(end - start) / CLOCKS_PER_SEC;

		// for each resolution in multi-resolution mesh:
		//   (1) print reading time, rendering time, (reading + rendering) time
		out.open(outinfo, ios::app);
		out << j << " \t" << reading_time << " \t" << rendering_time << " \t" << reading_time + rendering_time << endl;
		out.close();



		if (this->DisplayText)
			cout << "Level " << this->NumberOfFilters - j << ": " << this->SynthesisMeshes[j]->GetNumberOfCells()
			<< " faces, " << this->SynthesisMeshes[j]->GetNumberOfPoints() << " vertices and "
			<< this->SynthesisMeshes[j]->GetNumberOfEdges() << " edges " << endl;


//		if (this->Display != 0)
//		{
//			this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
//			this->MeshWindow->Render();
//
//#if ( (VTK_MAJOR_VERSION >= 5))
//			start2 = this->Timer->GetUniversalTime();
//			finish2 = this->Timer->GetUniversalTime();
//			while (finish2 - start2 < this->Time)
//				finish2 = this->Timer->GetUniversalTime();
//
//#else
//			start2 = this->Timer->GetCurrentTime();
//			finish2 = this->Timer->GetCurrentTime();
//			while (finish2 - start2 < this->Time)
//				finish2 = this->Timer->GetCurrentTime();
//#endif
//			if (this->Capture == 1)
//			{
//				std::stringstream strfile;
//				strfile << "Mesh" << this->NumberOfFilters - j << ".bmp";
//				this->MeshWindow->Capture(strfile.str().c_str());
//			}
//
//			if (this->Display == 2)
//			{
//				cout << "Window interaction: presse 'e' key to exit from interaction" << endl;
//				this->MeshWindow->Interact();
//			}
//		}
	}

	this->ArithmeticCoder->CloseFile();

	this->SynthesisMeshes[0]->Modified();

	// for each resolution in multi-resolution mesh:
	//   (2) save mesh as obj file

	// output file name initialization
	std::string strout;
	strout = this->InputFileName;
	if (this->Lifting > 0)
		strout += "_DS_lift_0_Mesh";
	else
		strout += "_DS_nolift_Mesh";

	double Factor, Tx, Ty, Tz;
	this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetScalingFactors(Factor, Tx, Ty, Tz);
	for (j = 0; j < this->NumberOfFilters; j++)
	{
		this->Filters[j]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
		this->Filters[j]->GetOutput()->UnQuantizeCoordinates();

		std::stringstream strfile;
		strfile << strout << this->NumberOfFilters - j << ".ply";
		vtkPLYWriter* Writer = vtkPLYWriter::New();
		Writer->SetFileName(strfile.str().c_str());
		Writer->SetInputData(this->Filters[j]->GetOutput());
		Writer->Write();
		Writer->Delete();
	}

	//if (this->Display != 0)
	//{
	//	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
	//	this->MeshWindow->Render();
	//	cout << "Window interaction: presse 'e' key to exit from interaction" << endl;
	//	this->MeshWindow->Interact();
	//}
	//if (this->WriteOutput == 1)
	//{
	//	// output file name initialization
	//	std::string strout;
	//	if (this->InputFileName != NULL)
	//	{
	//		strout = this->InputFileName;
	//		strout += "_";
	//	}
	//	else
	//		strout = "Mesh";

	//	this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->UnQuantizeCoordinates();
	//	if (this->FileType == 0)
	//	{
	//		std::stringstream strfile;
	//		strfile << strout << "0.iv";
	//		this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->WriteInventor(strfile.str().c_str());
	//	}
	//	else
	//	{
	//		vtkPLYWriter* Writer = vtkPLYWriter::New();
	//		std::stringstream strfile;
	//		strfile << strout << "0.ply";
	//		Writer->SetFileName(strfile.str().c_str());
	//		Writer->SetInputData(this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput());
	//		Writer->Write();
	//		Writer->Delete();
	//	}
	//	double Factor, Tx, Ty, Tz;
	//	this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetScalingFactors(Factor, Tx, Ty, Tz);
	//	for (j = 0; j < this->NumberOfFilters; j++)
	//	{
	//		this->Filters[j]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
	//		this->Filters[j]->GetOutput()->UnQuantizeCoordinates();
	//		std::stringstream strfile;
	//		if (this->FileType == 0)
	//		{
	//			strfile << strout << this->NumberOfFilters - j << ".iv";
	//			this->Filters[j]->GetOutput()->WriteInventor(strfile.str().c_str());
	//		}
	//		else
	//		{
	//			strfile << strout << this->NumberOfFilters - j << ".ply";
	//			vtkPLYWriter* Writer = vtkPLYWriter::New();
	//			Writer->SetFileName(strfile.str().c_str());
	//			Writer->SetInputData(this->Filters[j]->GetOutput());
	//			Writer->Write();
	//			Writer->Delete();
	//		}
	//	}
	//}
	//if (this->WriteOutput == 2)
	//{
	//	// output file name initialization
	//	std::string strout;
	//	if (this->InputFileName != NULL)
	//	{
	//		strout = this->InputFileName;
	//		strout += "_";
	//	}
	//	else
	//		strout = "Mesh";

	//	double Factor, Tx, Ty, Tz;
	//	this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetScalingFactors(Factor, Tx, Ty, Tz);
	//	if (this->FileType == 0)
	//	{
	//		this->Filters[0]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
	//		this->Filters[0]->GetOutput()->UnQuantizeCoordinates();

	//		std::stringstream strfile;
	//		strfile << strout << "0.iv";
	//		this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->WriteInventor(strfile.str().c_str());
	//	}
	//	else
	//	{
	//		for (j = 0; j < MIN(10, this->NumberOfFilters); j++)
	//		{
	//			this->Filters[j]->GetOutput()->SetScalingFactors(Factor, Tx, Ty, Tz);
	//			this->Filters[j]->GetOutput()->UnQuantizeCoordinates();

	//			std::stringstream strfile;
	//			if (this->FileType == 0)
	//			{
	//				strfile << strout << this->NumberOfFilters - j << ".iv";
	//				this->Filters[j]->GetOutput()->WriteInventor(strfile.str().c_str());
	//			}
	//			else
	//			{
	//				strfile << strout << this->NumberOfFilters - j << ".ply";
	//				vtkPLYWriter* Writer = vtkPLYWriter::New();
	//				Writer->SetFileName(strfile.str().c_str());
	//				Writer->SetInputData(this->Filters[j]->GetOutput());
	//				Writer->Write();
	//				Writer->Delete();
	//			}
	//		}
	//	}
	//}
}

void vtkMultiresolutionIO::Write()
{
	vtkIdType Id;
	vtkSurface *BaseMesh;

	this->PointsIds->SetNumberOfIds(this->Filters[0]->GetOutput()->GetNumberOfPoints());
	for (Id=0;Id<this->Filters[0]->GetOutput()->GetNumberOfPoints();Id++)
	{
		this->PointsIds->SetId(Id,Id);
	}
	BaseMesh=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput();

	this->ArithmeticCoder=vtkArithmeticCoder::New();
	this->ArithmeticCoder->OpenFile(this->FileName,1);
	this->ArithmeticCoder->InitConnectivityQSModels(1);
	this->ArithmeticCoder->StartCoding();
	this->ArithmeticCoder->EncodeBit(this->ArithmeticType!=0);
	this->EncodeMeshConnectivity(BaseMesh);
	this->EncodeMeshGeometry(BaseMesh);
	this->EncodeLiftingProperties();
	this->ArithmeticCoder->EncodeByte(this->NumberOfFilters);

	int Type,i;
	int bit1,bit2;
	for (i=this->NumberOfFilters-1;i>=0;i--)
	{
		Type=this->Filters[i]->GetSubdivisionType();
		bit1=Type&2;
		bit2=Type&1;
		this->ArithmeticCoder->EncodeBit(bit1!=0);
		this->ArithmeticCoder->EncodeBit(bit2!=0);
	}


	this->BaseMeshBitRate=this->ArithmeticCoder->StopCoding()*8;
	if (this->ArithmeticType==1)
		this->EncodeProgressiveResolution();
	else
		this->EncodeProgressivePrecision();
}
void vtkMultiresolutionIO::EncodeProgressiveResolution()
{
	int j;
	int GeometryBitrate[1000];
	int ConnectivityBitrate[1000];


	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		this->ArithmeticCoder->StartCoding();
		this->Filters[j]->SetIOType(1);
		this->Filters[j]->ArithmeticCoder=this->ArithmeticCoder;
		this->Filters[j]->SetPointsIds(this->PointsIds);
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);

		
		if (j==this->NumberOfFilters-1)
		{
			// *** We modify MergeInput's scaling factor... so this cannot be used
			//if (this->Filters[0]->GetMergeInput()!=0)
   //             this->EncodeScalingFactors(this->Filters[0]->GetMergeInput());
			//else
			//	this->EncodeScalingFactors(this->Filters[0]->GetOutput());
			this->EncodeScalingFactors(this->Filters[0]->GetOutput());
		}

		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();

		this->Filters[j]->Subdivide();
		ConnectivityBitrate[j]=this->ArithmeticCoder->StopCoding()*8;
		this->ArithmeticCoder->StartCoding();
		this->Filters[j]->WriteCoefficients();
		this->Filters[j]->Reconstruct();	
		GeometryBitrate[j]=this->ArithmeticCoder->StopCoding()*8;
	}
	this->ArithmeticCoder->CloseFile();

	if (this->WriteRepport > 0 && this->WriteRepport < 4)
	{
		int connectivitybits=0;
		double encodedbits=0,duration,ratio;

#if ( (VTK_MAJOR_VERSION >= 5))
		duration = this->Timer->GetUniversalTime()- this->StartTime;
#else
		duration = this->Timer->GetCurrentTime()- this->StartTime;
#endif
		ratio= (this->Filters[0]->GetOutput()->GetNumberOfCells())/duration;

		// for test2, test3
		std::string outinfo = this->InputFileName;
		if (this->WriteRepport == 2)
		{
			outinfo += "_lifting_repport.txt";
			std::ofstream out;
		}
		else if (this->WriteRepport == 3)
		{
			outinfo += "_lifting_random_repport.txt";
			std::ofstream out;
		}
		else
		{
			outinfo = "repport.txt";
		}
		std::ofstream Repport;
		Repport.open (outinfo, ofstream::out | ofstream::trunc);
		Repport<<"Filename: "<<this->FileName<<endl;

		Repport<<"Quantization : "<<this->Quantization<<" bits"<<endl;
		if (this->Lifting==0)
			Repport<<"No lifting"<<endl;
		else
		{
			if (this->Lifting==2)
				Repport<<"Fast 0-ring Lifting"<<endl;
			else
				Repport<<"Lifting radius: "<<this->LiftingRadius<<endl;
		}

		if (this->GeometricConstraint==0)
			Repport<<"No Wavelet Geometrical Criterion"<<endl;
		else
		{
			Repport<<"Geometry threshold: "<<this->EdgeAngleThreshold<<endl;
			Repport<<"Wavelet threshold: "<<this->WGC<<endl;
		}
		Repport<<"Total execution time : "<<duration<<" seconds : "<<ratio<<" faces/s"<<endl;


		double numberoffaces=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfCells();
		double numberofvertices=this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints();
		double cost=32.0+ceil(log(numberofvertices)/log(2.0))*3.0*numberoffaces;
		connectivitybits=(int) cost;
		encodedbits=this->BaseMeshBitRate;

		Repport<<"Level 0 : "<<this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfCells()
			<<"f, "<<this->Filters[this->NumberOfFilters-1]->GetSubdivisionInput()->GetNumberOfPoints()
			<<"v, total data: "<<(int) this->BaseMeshBitRate<<" bits (connectivity: "<<connectivitybits<<
			"bits)"<<endl;
		for (j=this->NumberOfFilters-1;j>=0;j--)
		{
			encodedbits+=ConnectivityBitrate[j]+GeometryBitrate[j];
			
			double RelativeBitrate=encodedbits/(double) this->Filters[0]->GetOutput()->GetNumberOfPoints();

			connectivitybits+=ConnectivityBitrate[j];
			Repport<<"Level "<<this->NumberOfFilters-j
				<<": "<<this->Filters[j]->GetOutput()->GetNumberOfCells()
				<<"f, "<<this->Filters[j]->GetOutput()->GetNumberOfPoints()
				<<"v, valence entropy= "<<this->Filters[j]->GetOutput()->GetValenceEntropy()
				<<", total data: "<<RelativeBitrate
				<<" bits/v (connectivity: "<<connectivitybits
				<<"bits, "<<
				(float)((float)ConnectivityBitrate[j])
				/((float)this->Filters[j]->GetOutput()->GetNumberOfPoints()-(float)this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints())
				<<" bits/vertex for this level)"<<endl;
		};

		Repport<<"Global coding: "<<(double) encodedbits/this->Filters[0]->GetOutput()->GetNumberOfPoints()
			<<" bits/vertex, connectivity : "<<
			(double) connectivitybits/this->Filters[0]->GetOutput()->GetNumberOfPoints()<<
			" bits/vertex, geometry : "<<
			(double) (encodedbits-connectivitybits)/this->Filters[0]->GetOutput()->GetNumberOfPoints()
			<<"bits/vertex"<<endl;
		Repport<<"File size: "<<((int) encodedbits)/8<<"bytes"<<endl;

		Repport.close();
	}

	// for DS project
	if (this->WriteRepport == 4)
	{
		int connectivitybits = 0;
		double encodedbits = 0, duration, ratio;

#if ( (VTK_MAJOR_VERSION >= 5))
		duration = this->Timer->GetUniversalTime() - this->StartTime;
#else
		duration = this->Timer->GetCurrentTime() - this->StartTime;
#endif
		ratio = (this->Filters[0]->GetOutput()->GetNumberOfCells()) / duration;

		std::string outinfo = this->InputFileName;
		if (this->Lifting > 0)
			outinfo += "_DS_lift_0_repport.txt";
		else
			outinfo += "_DS_nolift_repport.txt";
		std::ofstream Repport;
		Repport.open(outinfo, ofstream::out | ofstream::trunc);

		double numberoffaces = this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetNumberOfCells();
		double numberofvertices = this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetNumberOfPoints();
		double cost = 32.0 + ceil(log(numberofvertices) / log(2.0)) * 3.0 * numberoffaces;
		connectivitybits = (int)cost;
		encodedbits = this->BaseMeshBitRate;

		// level, v, f, base_connectivity, base_geometry, base_total
		Repport << "level \t" << "v \t" << "f \t" << "base_connectivity \t" << "base_geometry \t" << "base_total" << endl;
		Repport << this->NumberOfFilters
			<< " \t" << this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetNumberOfPoints()
			<< " \t" << this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetNumberOfCells()
			<< " \t" << connectivitybits
			<< " \t" << (int)this->BaseMeshBitRate - connectivitybits
			<< " \t" << (int)this->BaseMeshBitRate << endl;
		for (j = this->NumberOfFilters - 1; j >= 0; j--)
		{
			encodedbits += ConnectivityBitrate[j] + GeometryBitrate[j];
			connectivitybits += ConnectivityBitrate[j];

			// level, v, f, addt_connectivity, addt_geometry, addt_total
			Repport << j
				<< " \t" << this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints()
				<< " \t" << this->Filters[j]->GetSubdivisionInput()->GetNumberOfCells()
				<< " \t" << (int)ConnectivityBitrate[j]
				<< " \t" << (int)GeometryBitrate[j]
				<< " \t" << (int)ConnectivityBitrate[j] + (int)GeometryBitrate[j] << endl;
		};

		// global_connectivity, global_geometry, global_total
		Repport << (int)connectivitybits
			<< " \t" << (int)(encodedbits - connectivitybits)
			<< " \t" << encodedbits << endl;

		Repport.close();
	}
}

void vtkMultiresolutionIO::Analyse()
{
	double start2, finish2;
	int i,j,revert;

	this->AnalysisMeshes[0]=this->Input;
	this->MeshWindow->SetInputData(this->AnalysisMeshes[0]);
	if (this->Display==1)
	{
		this->MeshWindow->Render();
		cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
		this->MeshWindow->Interact();
	}

#if ( (VTK_MAJOR_VERSION >= 5))
	this->StartTime = this->Timer->GetUniversalTime();
#else
	this->StartTime = this->Timer->GetCurrentTime();
#endif

	i=0;
	int merged=1;
	while (merged==1 && i < 2)
	{
		std::vector<vtkIdType> bad_faces;
		bool badfacefull_flag = false;

		if (i+1>=this->MaxNumberOfLevels)
		{
			merged=0;
		}
		else
		{
			this->Filters[i]=vtkWaveletSubdivisionFilter::New();
			this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
			this->Filters[i]->SetLifting(this->Lifting);
			this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
			this->Filters[i]->SetArithmeticType(this->ArithmeticType);
			this->Filters[i]->Quantization=this->Quantization;
			this->Filters[i]->SetMergeInput(this->AnalysisMeshes[i]);
			this->Filters[i]->SetGeometryCriterion(this->GeometricConstraint);
			this->Filters[i]->GeometryPrediction=this->GeometryPrediction;
			this->Filters[i]->SetCurvatureTreshold(this->EdgeAngleThreshold);
			this->Filters[i]->SetWaveletTreshold(this->WGC);
			this->Filters[i]->SetDisplayEfficiency(this->DisplayEfficiency);

			this->Filters[i]->bad_faces = bad_faces;

			this->Filters[i]->SolveInverseProblem(0);
			this->AnalysisMeshes[i+1]=this->Filters[i]->GetMergeOutput();

			j = 0; // number of times of which bad faces become full

			while (((this->AnalysisMeshes[i+1]->GetNumberOfPoints()==this->AnalysisMeshes[i]->GetNumberOfPoints())
				||(this->Filters[i]->remaining_faces!=0))&&(j<MIN(10, this->AnalysisMeshes[i]->GetNumberOfCells() - 1)))//(j<this->AnalysisMeshes[i]->GetNumberOfCells()-1))
			{
				//cout << "Fail " << j << ", " << this->Filters[i]->remaining_faces << " remaining faces" << endl;
				bad_faces = this->Filters[i]->bad_faces;
				if (bad_faces.size() > 0)
				{
					/*if (bad_faces.size() == this->Filters[i]->GetMergeInput()->GetNumberOfCells())
					{
						bad_faces.clear();
						j++;
					}*/
					j--;
				}

				this->Filters[i]->Delete();
				this->Filters[i]=vtkWaveletSubdivisionFilter::New();
				this->Filters[i]->SetMergeInput(this->AnalysisMeshes[i]);
				this->Filters[i]->SetGeometryPrediction(this->GeometryPrediction);
				this->Filters[i]->SetLifting(this->Lifting);
				this->Filters[i]->SetLiftingRadius(this->LiftingRadius);
				this->Filters[i]->Quantization=this->Quantization;
				this->Filters[i]->SetArithmeticType(this->ArithmeticType);
				this->Filters[i]->SetGeometryCriterion(this->GeometricConstraint);
				this->Filters[i]->SetCurvatureTreshold(this->EdgeAngleThreshold);
				this->Filters[i]->SetWaveletTreshold(this->WGC);
				this->Filters[i]->SetDisplayEfficiency(this->DisplayEfficiency);

				this->Filters[i]->bad_faces = bad_faces;

				j++;
//				this->AnalysisMeshes[i+1]->Delete();
				this->Filters[i]->SolveInverseProblem(j);
				this->AnalysisMeshes[i+1]=this->Filters[i]->GetMergeOutput();
			}

			if (this->DisplayText && this->Filters[i]->remaining_faces == 0)
			{
				cout<<"i="<<i<<" ,original : "<<this->Filters[i]->GetMergeInput()->GetNumberOfCells()<<" faces, "<<this->Filters[i]->GetMergeInput()->GetPoints()->GetNumberOfPoints()<<" vertices";
				cout<<", simplified : "<<this->Filters[i]->GetMergeOutput()->GetNumberOfCells()<<" faces, "<<this->Filters[i]->GetMergeOutput()->GetPoints()->GetNumberOfPoints()<<" vertices, "<<this->Filters[i]->GetMergeOutput()->GetNumberOfEdges()<<" edges"<<endl;
			}

			if (this->Filters[i]->remaining_faces!=0)
				revert=0;

			if (this->Display==1)
			{
				if (this->GeometricConstraint==1)
				{	
					vtkLookupTable *lut=vtkLookupTable::New();
					this->MeshWindow->SetLookupTable(lut);
					this->AnalysisMeshes[i+1]->ComputeSharpVertices(this->EdgeAngleThreshold);
					lut->SetHueRange(0.667,0.0);
				}
				this->MeshWindow->SetInputData(this->AnalysisMeshes[i+1]);
				this->MeshWindow->Render();
#if ( (VTK_MAJOR_VERSION >= 5))
				start2=this->Timer->GetUniversalTime();
				finish2=this->Timer->GetUniversalTime();
				while (finish2 - start2<this->Time)
					finish2=this->Timer->GetUniversalTime();

#else
				start2=this->Timer->GetCurrentTime();
				finish2=this->Timer->GetCurrentTime();
				while (finish2 - start2<this->Time)
					finish2=this->Timer->GetCurrentTime();
#endif
				if (this->Display==1)
				{
					cout<<"Window interaction: presse 'e' key to exit from interaction"<<endl;
					this->MeshWindow->Interact();
				}
			}

			if ((this->AnalysisMeshes[i]->GetNumberOfCells()>0)
				&&(this->Filters[i]->GetMergeOutput()->GetNumberOfPoints()<this->Filters[i]->GetMergeInput()->GetNumberOfPoints()))
				merged=1;
			else
				merged=0;
			if (this->Filters[i]->remaining_faces!=0)
				merged=0;
		}
		i++;

	}
	this->NumberOfFilters=i-1;


	if (this->Display == 1)
	{
		this->MeshWindow->Interact();
	}
}
void vtkMultiresolutionIO::Synthetize()
{
	double start2, finish2;
	vtkIdType Id;
	int j;

	this->PointsIds->SetNumberOfIds(this->AnalysisMeshes[0]->GetNumberOfPoints());
	for (Id=0;Id<this->AnalysisMeshes[0]->GetNumberOfPoints();Id++)
	{
		this->PointsIds->SetId(Id,Id);
	}

	this->SynthesisMeshes[this->NumberOfFilters]=vtkSurface::New();
	this->SynthesisMeshes[this->NumberOfFilters]->CreateFromPolyData(this->AnalysisMeshes[this->NumberOfFilters]);

	double Factor, Tx, Ty, Tz;
	this->Filters[0]->GetMergeInput()->GetScalingFactors(Factor, Tx, Ty, Tz);
	this->SynthesisMeshes[this->NumberOfFilters]->SetScalingFactors(Factor, Tx, Ty, Tz);

	for (j=this->NumberOfFilters-1;j>=0;j--)
	{
		this->Filters[j]->SetPointsIds(this->PointsIds);
		this->Filters[j]->SetInput(this->SynthesisMeshes[j+1]);			
		this->Filters[j]->Subdivide();

		this->SynthesisMeshes[j]=this->Filters[j]->GetOutput();
		this->SynthesisMeshes[j]->SetScalingFactors(Factor, Tx, Ty, Tz);

		if (this->Display==1)
		{
			this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
			this->MeshWindow->Render();
			this->MeshWindow->Interact();
#if ( (VTK_MAJOR_VERSION >= 5))
			start2=this->Timer->GetUniversalTime();
			finish2=this->Timer->GetUniversalTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetUniversalTime();

#else
			start2=this->Timer->GetCurrentTime();
			finish2=this->Timer->GetCurrentTime();
			while (finish2 - start2<this->Time)
				finish2=this->Timer->GetCurrentTime();
#endif
		}
	}
	if (this->Display==1)
	{
		this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
		this->MeshWindow->Render();
		this->MeshWindow->Interact();
	}

	//// test : delete AnalysisMeshes[0~last-1] for reducing memory
	//for (j = this->NumberOfFilters - 1; j >= 0; j--)
	//{
	//	this->AnalysisMeshes[j]->Delete();
	//	this->Filters[j]->InitMergeInput();
	//}
}

void vtkMultiresolutionIO::DisplayHires()
{
	this->MeshWindow->SetInputData(this->SynthesisMeshes[this->NumberOfFilters-1]);
	this->MeshWindow->Render();
	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);
	this->MeshWindow->Render();
	this->MeshWindow->Interact();
}	

vtkMultiresolutionIO* vtkMultiresolutionIO::New()
{
	// First try to create the object from the vtkObjectFactory
	vtkObject* ret = vtkObjectFactory::CreateInstance("vtkMultiresolutionIO");
	if(ret)
	{
		return (vtkMultiresolutionIO*)ret;
	}
	// If the factory was unable to create the object, then create it here.
	return (new vtkMultiresolutionIO);
}

void vtkMultiresolutionIO::RoiSetting()
{
	int i,j;
	vtkIdType PtId;

	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);

	// Reloaded exported mesh
	std::string strfile = this->InputFileName;
	strfile += "_m.obj";

	vtkSurface *Mesh2 = vtkSurface::New();
	cout << "Load : " << strfile << endl;
	Mesh2->CreateFromFile(strfile.c_str());
	Mesh2->DisplayMeshProperties();

	vtkMultiresolutionIO *MIO2 = vtkMultiresolutionIO::New();
	MIO2->SetDisplay(this->Display);
	MIO2->SetArithmeticType(this->ArithmeticType);
	MIO2->SetQuantization(this->Quantization);
	MIO2->SetNumberOfBitPlanes(this->NumberOfBitPlanes);
	if (this->Lifting == 0)
	{
		MIO2->SetLifting(0);
	}
	else
	{
		MIO2->SetLifting(1);
		MIO2->SetLiftingRadius(this->LiftingRadius);
		cout << "Lifting radius : " << this->LiftingRadius << endl;
	}
	MIO2->SetGeometricalConstraint(this->GeometricConstraint);
	MIO2->SetEdgeAngleThreshold(this->EdgeAngleThreshold);
	MIO2->SetWGC(this->WGC);

	if (this->ArithmeticType == 1)
		Mesh2->QuantizeCoordinatesLike(this->Filters[0]->GetOutput());

	MIO2->SetInput(Mesh2);
	MIO2->SetInput2();

	// sieve for modified vertices (init with all 0s)
	bool *ModifiedTable = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
	int number_of_modified = 0;

	////////// compare and replace orig with modified ones //////////
	{
		vtkSurface *Surface1 = this->Filters[0]->GetOutput();
		vtkSurface *Surface2 = (vtkSurface *)(MIO2->GetRenderWindow()->GetInput());

		///// # of v and # of f  equality test /////
		if ((Surface1->GetNumberOfPoints() != Surface2->GetNumberOfPoints())
			|| (Surface1->GetNumberOfCells() != Surface2->GetNumberOfCells()))
		{
			// for debug //
			cout << "Compare: " << endl;
			cout << "mesh1 has " << Surface1->GetNumberOfPoints() << " vertices, " << Surface1->GetNumberOfCells() << " faces" << endl;
			cout << "mesh2 has " << Surface2->GetNumberOfPoints() << " vertices, " << Surface2->GetNumberOfCells() << " faces" << endl;
		}
		else
		{
			vtkPoints *Points1 = Surface1->GetPoints();
			vtkPoints *Points2 = Surface2->GetPoints();

			vtkIdType i;
			double e = 1e-3;
			vtkPoints *BasePoints = this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetPoints();
			int BaseSize = this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetNumberOfPoints();

			for (i = 0; i < Surface1->GetNumberOfPoints(); i++)
			{
				double p1[3], p2[3];
				Points1->GetPoint(i, p1);
				Points2->GetPoint(i, p2);
				if (abs(p1[0] - p2[0]) > e || abs(p1[1] - p2[1]) > e || abs(p1[2] - p2[2]) > e)
				{
					number_of_modified++;

					// Replace with new coords and add ROI
					// 
					// WARNING!! : Mesh coords are not propagated to other meshes (=AnalysisMeshes[],SynthesisMeshes[])
					//   The coords required to be written to file are 'Basemesh' and 'Wavelets'.
					//   Basemesh is stored in SynthesisMeshes[last] so we need to change it also.

					if (i < BaseSize)
						BasePoints->SetPoint(i, p2);  // for storing basemesh to file

					Points1->SetPoint(i, p2);
					ModifiedTable[i] = true;
				}
			}
			BasePoints->Modified();
			Points1->Modified();
		}
	}
	Mesh2->Delete();
	MIO2->Delete();

	cout << "Start ROI propagation" << endl;

	clock_t start, end;

	// ROI propagation
	// 1) without lifting
	//		ROI must be propagated since w = v - (p1+p2)/2 and v,p1,p2 were modified
	//		v is recalculated from user input.
	//		if v is not remained at last level, v turns into w, so w_old must be recalculated
	//		if v is remained at last level, no change required
	//    Optimization:
	//      for each wavelet, check if v/p1/p2 have been changed and recalculate wavelet.
	if (this->Lifting == 0)
	{
		// for test2, test3
		std::string outinfo = this->InputFileName;
		outinfo = outinfo.substr(0, outinfo.size() - 4);
		outinfo += "_info.txt";
		std::ofstream out;

		start = clock();

		double V[3], P1[3], P2[3];
		int number_of_handled = 0;
		vtkIdType PtId;
		vtkPoints *Points = this->Filters[0]->GetOutput()->GetPoints();

		// Roi wavelets modification
		if (this->ArithmeticType != 0)
		{
			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
				vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
				int *IntWavelet;

				// current range propagate
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					vtkIdType p1 = this->Filters[j]->vertices[PtId].parent1;
					vtkIdType p2 = this->Filters[j]->vertices[PtId].parent2;

					if (ModifiedTable[PtId] == true || ModifiedTable[p1] == true || ModifiedTable[p2] == true)
					{
						number_of_handled++;
						Points->GetPoint(PtId, V);
						Points->GetPoint(p1, P1);
						Points->GetPoint(p2, P2);
						IntWavelet = this->Filters[j]->IntegerWavelets->GetPointer((PtId - NOldPoints) * 3);

						IntWavelet[0] = (int)floor(V[0] - 0.5*(P1[0] + P2[0]));
						IntWavelet[1] = (int)floor(V[1] - 0.5*(P1[1] + P2[1]));
						IntWavelet[2] = (int)floor(V[2] - 0.5*(P1[2] + P2[2]));
					}
				}
				this->Filters[j]->IntegerWavelets->Modified();
			}
		}
		else
		{
			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
				vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
				double *Wavelet;

				// current range propagate
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					vtkIdType p1 = this->Filters[j]->vertices[PtId].parent1;
					vtkIdType p2 = this->Filters[j]->vertices[PtId].parent2;

					if (ModifiedTable[PtId] == true || ModifiedTable[p1] == true || ModifiedTable[p2] == true)
					{
						number_of_handled++;
						Points->GetPoint(PtId, V);
						Points->GetPoint(p1, P1);
						Points->GetPoint(p2, P2);
						Wavelet = this->Filters[j]->Wavelets->GetPointer((PtId - NOldPoints) * 3);

						Wavelet[0] = V[0] - 0.5*(P1[0] + P2[0]);
						Wavelet[1] = V[1] - 0.5*(P1[1] + P2[1]);
						Wavelet[2] = V[2] - 0.5*(P1[2] + P2[2]);
					}
				}
				this->Filters[j]->Wavelets->Modified();
			}
		}

		end = clock();

		// Basemesh (for calculate number of modified vertices in basemesh)
		{
			vtkIdType NNewPoints = this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput()->GetNumberOfPoints();
			vtkIdType NOldPoints = 0;
			for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
			{
				if (ModifiedTable[PtId] == true)
					number_of_handled++;
			}
		}

		// for test2
		out.open(outinfo, ios::app);
		out << "ROI wavelet modification elapsed time: " << double(end - start) / CLOCKS_PER_SEC << endl;
		out << "Number of Modified vertices : " << number_of_modified << endl;
		out << "Number of handled vertices: " << number_of_handled << endl;
		out.close();

		start = clock();
		// Newly quantizing
		//		No-lifting case(only...?) need vertex quantization and wavelet quantization
		//      w = floor((w_orig / factor_orig) + T_orig - T_new) * factor_new + 0.5
		//        = floor( w_orig * (factor_new / factor_orig) + (T_orig - T_new) * factor_new + 0.5 )
		if (this->ArithmeticType == 1)
		{
			double factor_orig, factor_new, T_orig[3], T_new[3];
			this->Filters[0]->GetOutput()->GetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);

			// unquantize vertices
			vtkSurface *OrigMesh = this->Filters[0]->GetOutput();
			vtkSurface *BaseMesh = this->Filters[this->NumberOfFilters - 1]->GetSubdivisionInput();
			OrigMesh->SetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);
			BaseMesh->SetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);
			OrigMesh->UnQuantizeCoordinates();
			BaseMesh->UnQuantizeCoordinates();
			OrigMesh->Modified();
			BaseMesh->Modified();

			// requantize vertices
			OrigMesh->QuantizeCoordinates(this->Quantization);
			OrigMesh->GetScalingFactors(factor_new, T_new[0], T_new[1], T_new[2]);
			BaseMesh->QuantizeCoordinatesLike(OrigMesh);
			OrigMesh->Modified();
			BaseMesh->Modified();

			// un-requantize wavelets
			double Wavelet[3];
			for (j = 0; j < this->NumberOfFilters - 1; j++)
			{
				vtkIntArray *IntWavelets = this->Filters[j]->IntegerWavelets;
				for (i = 0; i < IntWavelets->GetNumberOfTuples(); i++)
				{
					// IntWavelet = IntWavelets->GetPointer((NonZero->j) * 3);    is this maybe suitable?
					IntWavelets->GetTuple(i, Wavelet);
					Wavelet[0] = floor(Wavelet[0] * (factor_new / factor_orig) + (T_new[0] - T_orig[0]) * factor_orig + 0.5);
					Wavelet[1] = floor(Wavelet[1] * (factor_new / factor_orig) + (T_new[1] - T_orig[1]) * factor_orig + 0.5);
					Wavelet[2] = floor(Wavelet[2] * (factor_new / factor_orig) + (T_new[2] - T_orig[2]) * factor_orig + 0.5);
				}
				IntWavelets->Modified();
			}
		}
		else
		{
			;
		}
		end = clock();

		// for test2
		out.open(outinfo, ios::app);
		out << "Re-quantization elapsed time: " << double(end - start) / CLOCKS_PER_SEC << endl;
		out << "Number of Modified vertices : " << number_of_modified << endl;
		out.close();
	}
	
	// 2) with lifting
	//      after approximate,
	//		v is recalculated from Mesh2
	//		all vertices must be recalculated sequentially from finest to coarsest
	//      all SynthesisMeshes must be handled.
	if (this->Lifting == 1)
	{
		// for test2, test3
		std::string outinfo = this->FileName;
		outinfo = outinfo.substr(0, outinfo.size() - 7);
		outinfo += "info.txt";
		std::ofstream out;

		start = clock();
		// Newly quantizing
		//		lifting case also need vertex quantization
		//      all SynthesisMeshes must be requantized with same factor as MergeInput
		if (this->ArithmeticType == 1)
		{
			double factor_orig, factor_new, T_orig[3], T_new[3];
			this->Filters[0]->GetOutput()->GetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);

			vtkSurface *OrigMesh = this->Filters[0]->GetOutput();
			OrigMesh->SetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);
			OrigMesh->UnQuantizeCoordinates();
			OrigMesh->Modified();

			OrigMesh->QuantizeCoordinates(this->Quantization);
			OrigMesh->GetScalingFactors(factor_new, T_new[0], T_new[1], T_new[2]);
			OrigMesh->Modified();

			// un-re-quantize every mesh
			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkSurface *LowerMesh = this->Filters[j]->GetSubdivisionInput();
				LowerMesh->SetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);
				LowerMesh->UnQuantizeCoordinates();
				LowerMesh->Modified();

				LowerMesh->QuantizeCoordinates(factor_new, T_new[0], T_new[1], T_new[2]);
				LowerMesh->Modified();
			}
		}
		else
		{
			;
		}
		end = clock();

		// for test3
		out.open(outinfo, ios::app);
		out << "Re-quantization elapsed time: " << double(end - start) / CLOCKS_PER_SEC << endl;
		out.close();

		start = clock();

		// R^finest initialize
		bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
		for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
			R_raw[PtId] = ModifiedTable[PtId];

		vtkIdList *Edges = vtkIdList::New();
		vtkIdType PtId;
		vtkIdType edge;
		vtkIdType c, p1, p2;

		vtkIdList *HorizontalList, *VerticalList, *TempList;
		HorizontalList = vtkIdList::New();
		VerticalList = vtkIdList::New();

		vtkSparseMatrix::NonZeroElement *NonZero;

		for (j = 0; j < this->NumberOfFilters; j++)
		{
			vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
			vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
			vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
			vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

			bool *R_next = new bool[NOldPoints]();
			bool *W_next = new bool[NNewPoints-NOldPoints]();

			// w impacted by v in R^(j+1) must be added to hat(W)^j
			//     and propagated into k-disk of its parents
			for (PtId = 0; PtId < NOldPoints; PtId++)
			{
				if (R_raw[PtId] == true)
				{
					// add v to R^j 
					R_next[PtId] = true;

					this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
					for (i = 0; i < Edges->GetNumberOfIds(); i++)
					{
						// find c impacted by v
						edge = Edges->GetId(i);
						c = this->Filters[j]->edgesvector[edge].child;
						if ( (c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
						{
							// add c to hat(W)^j
							W_next[c - NOldPoints] = true;

							// propagate parents of c into k-disk
							HorizontalList->SetNumberOfIds(2);
							HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
							HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

							if (this->LiftingRadius > 0)
							{
								for (i = 0; i < this->LiftingRadius; i++)
								{
									this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
									TempList = VerticalList;
									VerticalList = HorizontalList;
									HorizontalList = TempList;
								}
							}

							// all v in k-disk are added to R^j
							int Range = HorizontalList->GetNumberOfIds();
							for (i = 0; i < Range; i++)
							{
								c = HorizontalList->GetId(i);
								if (c < NOldPoints)
									R_next[HorizontalList->GetId(i)] = true;
							}
						}
					}
				}
			}
			// mergable v in R^(j+1) must be added to hat(W)^j
			for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
			{
				if (R_raw[PtId] == true)
				{
					W_next[PtId - NOldPoints] = true;

					// propagate parents of PtId into k-disk
					HorizontalList->SetNumberOfIds(2);
					HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
					HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

					if (this->LiftingRadius > 0)
					{
						for (i = 0; i < this->LiftingRadius; i++)
						{
							this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
							TempList = VerticalList;
							VerticalList = HorizontalList;
							HorizontalList = TempList;
						}
					}

					// all v in k-disk are added to R^j
					int Range = HorizontalList->GetNumberOfIds();
					for (i = 0; i < Range; i++)
					{
						c = HorizontalList->GetId(i);
						if (c < NOldPoints)
							R_next[HorizontalList->GetId(i)] = true;
					}
				}
			}

			if (this->ArithmeticType == 1)
			{
				// recalculate w in hat(W)^j
				vtkIntArray *IntWavelets = this->Filters[j]->IntegerWavelets;
				int *IntWavelet;
				double P1[3], P2[3], V1[3], V[3];
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (W_next[PtId - NOldPoints] == true)
					{
						p1 = this->Filters[j]->vertices[PtId].parent1;
						p2 = this->Filters[j]->vertices[PtId].parent2;
						NewPoints->GetPoint(PtId, V1);
						NewPoints->GetPoint(p1, P1);
						NewPoints->GetPoint(p2, P2);

						IntWavelet = IntWavelets->GetPointer((PtId - NOldPoints) * 3);

						IntWavelet[0] = (int)floor(V1[0] - 0.5*(P1[0] + P2[0]));
						IntWavelet[1] = (int)floor(V1[1] - 0.5*(P1[1] + P2[1]));
						IntWavelet[2] = (int)floor(V1[2] - 0.5*(P1[2] + P2[2]));
					}
				}
				IntWavelets->Modified();

				// recalulate v in R^j
				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_next[PtId] == true)
					{
						V[0] = 0;
						V[1] = 0;
						V[2] = 0;

						NonZero = this->Filters[j]->Alpha->FirstHor[PtId];

						while (NonZero)
						{
							IntWavelet = IntWavelets->GetPointer((NonZero->j) * 3);

							V[0] += NonZero->Value*IntWavelet[0];
							V[1] += NonZero->Value*IntWavelet[1];
							V[2] += NonZero->Value*IntWavelet[2];

							NonZero = NonZero->NextHor;
						}

						NewPoints->GetPoint(PtId, P1);

						P2[0] = P1[0] + floor(V[0] + 0.5);
						P2[1] = P1[1] + floor(V[1] + 0.5);
						P2[2] = P1[2] + floor(V[2] + 0.5);

						OldPoints->SetPoint(PtId, P2);
					}
				}
			}

			R_raw = R_next;
		}
		end = clock();

		// for test3
		out.open(outinfo, ios::app);
		out << "ROI vertices modification elapsed time: " << double(end - start) / CLOCKS_PER_SEC << endl;
		out.close();


		///////////////////////////////////////////////////////////
		////////// count modified vertices and wavelets ///////////
		///////////////////////////////////////////////////////////
		{
			// number of vertices and wavelets to be modified
			int number_of_modified_wavelets = 0;
			int number_of_modified_vertices = 0;

			// R^finest initialize
			bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
			for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
			{
				R_raw[PtId] = ModifiedTable[PtId];
				if (R_raw[PtId] == true)
					number_of_modified_vertices += 1;
			}

			vtkIdList *Edges = vtkIdList::New();
			vtkIdType PtId;
			vtkIdType edge;
			vtkIdType c, p1, p2;

			vtkIdList *HorizontalList, *VerticalList, *TempList;
			HorizontalList = vtkIdList::New();
			VerticalList = vtkIdList::New();

			vtkSparseMatrix::NonZeroElement *NonZero;

			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
				vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
				vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
				vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

				bool *R_next = new bool[NOldPoints]();
				bool *W_next = new bool[NNewPoints - NOldPoints]();

				// w impacted by v in R^(j+1) must be added to hat(W)^j
				//     and propagated into k-disk of its parents
				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						// add v to R^j 
						R_next[PtId] = true;

						this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
						for (i = 0; i < Edges->GetNumberOfIds(); i++)
						{
							// find c impacted by v
							edge = Edges->GetId(i);
							c = this->Filters[j]->edgesvector[edge].child;
							if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
							{
								// add c to hat(W)^j
								W_next[c - NOldPoints] = true;

								// propagate parents of c into k-disk
								HorizontalList->SetNumberOfIds(2);
								HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
								HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

								if (this->LiftingRadius > 0)
								{
									for (i = 0; i < this->LiftingRadius; i++)
									{
										this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
										TempList = VerticalList;
										VerticalList = HorizontalList;
										HorizontalList = TempList;
									}
								}

								// all v in k-disk are added to R^j
								int Range = HorizontalList->GetNumberOfIds();
								for (i = 0; i < Range; i++)
								{
									c = HorizontalList->GetId(i);
									if (c < NOldPoints)
									{
										R_next[HorizontalList->GetId(i)] = true;
									}
								}
							}
						}
					}
				}
				// mergable v in R^(j+1) must be added to hat(W)^j
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						W_next[PtId - NOldPoints] = true;

						// propagate parents of PtId into k-disk
						HorizontalList->SetNumberOfIds(2);
						HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
						HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

						if (this->LiftingRadius > 0)
						{
							for (i = 0; i < this->LiftingRadius; i++)
							{
								this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
								TempList = VerticalList;
								VerticalList = HorizontalList;
								HorizontalList = TempList;
							}
						}

						// all v in k-disk are added to R^j
						int Range = HorizontalList->GetNumberOfIds();
						for (i = 0; i < Range; i++)
						{
							c = HorizontalList->GetId(i);
							if (c < NOldPoints)
							{
								R_next[HorizontalList->GetId(i)] = true;
							}
						}
					}
				}

				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_next[PtId] == true)
						number_of_modified_vertices += 1;
				}
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (W_next[PtId - NOldPoints] == true)
						number_of_modified_wavelets += 1;
				}

				R_raw = R_next;
			}

			int number_of_vertices = this->Filters[0]->GetOutput()->GetNumberOfPoints();

			out.open(outinfo, ios::app);
			out << " " << number_of_modified;
			out << " / " << number_of_vertices;
			out << " = " << (double)number_of_modified / number_of_vertices << "%" << endl;
			out << " " << number_of_modified_vertices;
			out << " " << number_of_modified_wavelets;
			out << " " << number_of_modified_vertices + number_of_modified_wavelets << endl;
			out.close();
		}
	}

	cout << "End ROI propagation" << endl;
}

// randomly modify n/k vertices
void vtkMultiresolutionIO::RoiRandom(int k)
{
	int i, j, e;
	vtkIdType PtId;

	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);

	// 2) with lifting
	//      after approximate,
	//		v is recalculated from Mesh2
	//		all vertices must be recalculated sequentially from finest to coarsest
	//      all SynthesisMeshes must be handled.
	if (this->Lifting == 1)
	{
		// propagation //
		{
			
			// Newly quantizing
			//		lifting case also need vertex quantization
			//      all SynthesisMeshes must be requantized with same factor as MergeInput
			if (this->ArithmeticType == 1)
			{
				double factor_orig, factor_new, T_orig[3], T_new[3];
				this->Filters[0]->GetOutput()->GetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);

				vtkSurface *OrigMesh = this->Filters[0]->GetOutput();
				OrigMesh->SetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);
				OrigMesh->UnQuantizeCoordinates();
				OrigMesh->Modified();

				OrigMesh->QuantizeCoordinates(this->Quantization);
				OrigMesh->GetScalingFactors(factor_new, T_new[0], T_new[1], T_new[2]);
				OrigMesh->Modified();

				// un-re-quantize every mesh
				for (j = 0; j < this->NumberOfFilters; j++)
				{
					vtkSurface *LowerMesh = this->Filters[j]->GetSubdivisionInput();
					LowerMesh->SetScalingFactors(factor_orig, T_orig[0], T_orig[1], T_orig[2]);
					LowerMesh->UnQuantizeCoordinates();
					LowerMesh->Modified();

					LowerMesh->QuantizeCoordinates(factor_new, T_new[0], T_new[1], T_new[2]);
					LowerMesh->Modified();
				}
			}
			else
			{
				;
			}

			// R^finest initialize
			bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
			for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
				R_raw[PtId] = true;

			vtkIdList *Edges = vtkIdList::New();
			vtkIdType PtId;
			vtkIdType edge;
			vtkIdType c, p1, p2;

			vtkIdList *HorizontalList, *VerticalList, *TempList;
			HorizontalList = vtkIdList::New();
			VerticalList = vtkIdList::New();

			vtkSparseMatrix::NonZeroElement *NonZero;

			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
				vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
				vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
				vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

				bool *R_next = new bool[NOldPoints]();
				bool *W_next = new bool[NNewPoints - NOldPoints]();

				int lengh_of_R_curr = 0;

				// w impacted by v in R^(j+1) must be added to hat(W)^j
				//     and propagated into k-disk of its parents
				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						// add v to R^j 
						R_next[PtId] = true;

						this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
						for (e = 0; e < Edges->GetNumberOfIds(); e++)
						{
							// find c impacted by v
							edge = Edges->GetId(e);
							c = this->Filters[j]->edgesvector[edge].child;
							if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
							{
								// add c to hat(W)^j
								W_next[c - NOldPoints] = true; //////////////////////////////////////// cannot enter this???

								// propagate parents of c into k-disk
								HorizontalList->SetNumberOfIds(2);
								HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
								HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

								if (this->LiftingRadius > 0)
								{
									for (i = 0; i < this->LiftingRadius; i++)
									{
										this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
										TempList = VerticalList;
										VerticalList = HorizontalList;
										HorizontalList = TempList;
									}
								}

								// all v in k-disk are added to R^j
								int Range = HorizontalList->GetNumberOfIds();
								for (i = 0; i < Range; i++)
								{
									c = HorizontalList->GetId(i);
									if (c < NOldPoints)
										R_next[HorizontalList->GetId(i)] = true;
								}
							}
						}
					}
				}
				// mergable v in R^(j+1) must be added to hat(W)^j
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						W_next[PtId - NOldPoints] = true;

						// propagate parents of PtId into k-disk
						HorizontalList->SetNumberOfIds(2);
						HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
						HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

						if (this->LiftingRadius > 0)
						{
							for (i = 0; i < this->LiftingRadius; i++)
							{
								this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
								TempList = VerticalList;
								VerticalList = HorizontalList;
								HorizontalList = TempList;
							}
						}

						// all v in k-disk are added to R^j
						int Range = HorizontalList->GetNumberOfIds();
						for (i = 0; i < Range; i++)
						{
							c = HorizontalList->GetId(i);
							if (c < NOldPoints)
								R_next[HorizontalList->GetId(i)] = true;
						}
					}
				}

				if (this->ArithmeticType == 1)
				{
					// recalculate w in hat(W)^j
					vtkIntArray *IntWavelets = this->Filters[j]->IntegerWavelets;
					int *IntWavelet;
					double P1[3], P2[3], V1[3], V[3];
					for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
					{
						if (W_next[PtId - NOldPoints] == true)
						{
							p1 = this->Filters[j]->vertices[PtId].parent1;
							p2 = this->Filters[j]->vertices[PtId].parent2;
							NewPoints->GetPoint(PtId, V1);
							NewPoints->GetPoint(p1, P1);
							NewPoints->GetPoint(p2, P2);

							IntWavelet = IntWavelets->GetPointer((PtId - NOldPoints) * 3);

							IntWavelet[0] = (int)floor(V1[0] - 0.5*(P1[0] + P2[0]));
							IntWavelet[1] = (int)floor(V1[1] - 0.5*(P1[1] + P2[1]));
							IntWavelet[2] = (int)floor(V1[2] - 0.5*(P1[2] + P2[2]));
						}
					}
					IntWavelets->Modified();

					// recalulate v in R^j
					for (PtId = 0; PtId < NOldPoints; PtId++)
					{
						if (R_next[PtId] == true)
						{
							V[0] = 0;
							V[1] = 0;
							V[2] = 0;

							NonZero = this->Filters[j]->Alpha->FirstHor[PtId];

							while (NonZero)
							{
								IntWavelet = IntWavelets->GetPointer((NonZero->j) * 3);

								V[0] += NonZero->Value*IntWavelet[0];
								V[1] += NonZero->Value*IntWavelet[1];
								V[2] += NonZero->Value*IntWavelet[2];

								NonZero = NonZero->NextHor;
							}

							NewPoints->GetPoint(PtId, P1);

							P2[0] = P1[0] + floor(V[0] + 0.5);
							P2[1] = P1[1] + floor(V[1] + 0.5);
							P2[2] = P1[2] + floor(V[2] + 0.5);

							OldPoints->SetPoint(PtId, P2);
						}
					}
				}

				R_raw = R_next;
			}
		}

		///////////////////////////////////////////////////////////
		////////// count modified vertices and wavelets ///////////
		///////////////////////////////////////////////////////////
		//{
		//	// number of vertices and wavelets to be modified
		//	int number_of_modified_wavelets = 0;
		//	int number_of_modified_vertices = 0;

		//	// R^finest initialize
		//	bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
		//	for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
		//	{
		//		R_raw[PtId] = ModifiedTable[PtId];
		//		if (R_raw[PtId] == true)
		//			number_of_modified_vertices += 1;
		//	}

		//	vtkIdList *Edges = vtkIdList::New();
		//	vtkIdType PtId;
		//	vtkIdType edge;
		//	vtkIdType c, p1, p2;

		//	vtkIdList *HorizontalList, *VerticalList, *TempList;
		//	HorizontalList = vtkIdList::New();
		//	VerticalList = vtkIdList::New();

		//	vtkSparseMatrix::NonZeroElement *NonZero;

		//	for (j = 0; j < this->NumberOfFilters; j++)
		//	{
		//		vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
		//		vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
		//		vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
		//		vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

		//		bool *R_next = new bool[NOldPoints]();
		//		bool *W_next = new bool[NNewPoints - NOldPoints]();

		//		// w impacted by v in R^(j+1) must be added to hat(W)^j
		//		//     and propagated into k-disk of its parents
		//		for (PtId = 0; PtId < NOldPoints; PtId++)
		//		{
		//			if (R_raw[PtId] == true)
		//			{
		//				// add v to R^j 
		//				R_next[PtId] = true;

		//				this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
		//				for (e = 0; e < Edges->GetNumberOfIds(); e++)
		//				{
		//					// find c impacted by v
		//					edge = Edges->GetId(e);
		//					c = this->Filters[j]->edgesvector[edge].child;
		//					if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
		//					{
		//						// add c to hat(W)^j
		//						W_next[c - NOldPoints] = true;

		//						// propagate parents of c into k-disk
		//						HorizontalList->SetNumberOfIds(2);
		//						HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
		//						HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

		//						if (this->LiftingRadius > 0)
		//						{
		//							for (i = 0; i < this->LiftingRadius; i++)
		//							{
		//								this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
		//								TempList = VerticalList;
		//								VerticalList = HorizontalList;
		//								HorizontalList = TempList;
		//							}
		//						}

		//						// all v in k-disk are added to R^j
		//						int Range = HorizontalList->GetNumberOfIds();
		//						for (i = 0; i < Range; i++)
		//						{
		//							c = HorizontalList->GetId(i);
		//							if (c < NOldPoints)
		//							{
		//								R_next[HorizontalList->GetId(i)] = true;
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}
		//		// mergable v in R^(j+1) must be added to hat(W)^j
		//		for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
		//		{
		//			if (R_raw[PtId] == true)
		//			{
		//				W_next[PtId - NOldPoints] = true;

		//				// propagate parents of PtId into k-disk
		//				HorizontalList->SetNumberOfIds(2);
		//				HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
		//				HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

		//				if (this->LiftingRadius > 0)
		//				{
		//					for (i = 0; i < this->LiftingRadius; i++)
		//					{
		//						this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
		//						TempList = VerticalList;
		//						VerticalList = HorizontalList;
		//						HorizontalList = TempList;
		//					}
		//				}

		//				// all v in k-disk are added to R^j
		//				int Range = HorizontalList->GetNumberOfIds();
		//				for (i = 0; i < Range; i++)
		//				{
		//					c = HorizontalList->GetId(i);
		//					if (c < NOldPoints)
		//					{
		//						R_next[HorizontalList->GetId(i)] = true;
		//					}
		//				}
		//			}
		//		}

		//		for (PtId = 0; PtId < NOldPoints; PtId++)
		//		{
		//			if (R_next[PtId] == true)
		//				number_of_modified_vertices += 1;
		//		}
		//		for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
		//		{
		//			if (W_next[PtId - NOldPoints] == true)
		//				number_of_modified_wavelets += 1;
		//		}

		//		R_raw = R_next;
		//	}

		//	cout << " " << (this->Filters[0]->GetOutput()->GetNumberOfPoints()) / k;
		//	cout << " " << number_of_modified_vertices + number_of_modified_wavelets << endl;
		//}
	}

	//cout << "End ROI propagation : " << (this->Filters[0]->GetOutput()->GetNumberOfPoints()) / k << endl;
}

// randomly modify n/k, 2n/k, ..., n vertices
void vtkMultiresolutionIO::RoiIncrease(int k)
{
	int i, j, t, e;
	vtkIdType PtId;

	std::string outinfo = this->InputFileName;
	outinfo += "_info.txt";
	std::ofstream out;

	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);

	// random number generator without duplicating
	std::vector<int> RandomNumbers(this->Filters[0]->GetOutput()->GetNumberOfPoints());
	for (int i = 0; i < RandomNumbers.size(); i++)
		RandomNumbers[i] = i;
	std::shuffle(RandomNumbers.begin(), RandomNumbers.end(), std::default_random_engine(std::random_device()()));

	if (k > this->Filters[0]->GetOutput()->GetNumberOfPoints())
		k = this->Filters[0]->GetOutput()->GetNumberOfPoints();

	int step = this->Filters[0]->GetOutput()->GetNumberOfPoints() / k;
	for (t = step; t < this->Filters[0]->GetOutput()->GetNumberOfPoints(); t += step)
	{
		// sieve for modified vertices (init with all 0s)
		bool *ModifiedTable = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();

		// randomly choose t vertices
		for (i = 0; i < t; i++)
			ModifiedTable[RandomNumbers[i]] = true;

		clock_t start, end;

		// 2) with lifting
		//      after approximate,
		//		v is recalculated from Mesh2
		//		all vertices must be recalculated sequentially from finest to coarsest
		//      all SynthesisMeshes must be handled.
		if (this->Lifting == 1)
		{
			// propagation //
			{
				start = clock();

				// R^finest initialize
				bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
				for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
					R_raw[PtId] = ModifiedTable[PtId];

				vtkIdList *Edges = vtkIdList::New();
				vtkIdType PtId;
				vtkIdType edge;
				vtkIdType c, p1, p2;

				vtkIdList *HorizontalList, *VerticalList, *TempList;
				HorizontalList = vtkIdList::New();
				VerticalList = vtkIdList::New();

				vtkSparseMatrix::NonZeroElement *NonZero;

				for (j = 0; j < this->NumberOfFilters; j++)
				{
					vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
					vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
					vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
					vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

					bool *R_next = new bool[NOldPoints]();
					bool *W_next = new bool[NNewPoints - NOldPoints]();

					int lengh_of_R_curr = 0;

					// w impacted by v in R^(j+1) must be added to hat(W)^j
					//     and propagated into k-disk of its parents
					for (PtId = 0; PtId < NOldPoints; PtId++)
					{
						if (R_raw[PtId] == true)
						{
							// add v to R^j 
							R_next[PtId] = true;

							this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
							for (e = 0; e < Edges->GetNumberOfIds(); e++)
							{
								// find c impacted by v
								edge = Edges->GetId(e);
								c = this->Filters[j]->edgesvector[edge].child;
								if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
								{
									// add c to hat(W)^j
									W_next[c - NOldPoints] = true; //////////////////////////////////////// cannot enter this???

									// propagate parents of c into k-disk
									HorizontalList->SetNumberOfIds(2);
									HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
									HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

									if (this->LiftingRadius > 0)
									{
										for (i = 0; i < this->LiftingRadius; i++)
										{
											this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
											TempList = VerticalList;
											VerticalList = HorizontalList;
											HorizontalList = TempList;
										}
									}

									// all v in k-disk are added to R^j
									int Range = HorizontalList->GetNumberOfIds();
									for (i = 0; i < Range; i++)
									{
										c = HorizontalList->GetId(i);
										if (c < NOldPoints)
											R_next[HorizontalList->GetId(i)] = true;
									}
								}
							}
						}
					}
					// mergable v in R^(j+1) must be added to hat(W)^j
					for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
					{
						if (R_raw[PtId] == true)
						{
							W_next[PtId - NOldPoints] = true;

							// propagate parents of PtId into k-disk
							HorizontalList->SetNumberOfIds(2);
							HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
							HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

							if (this->LiftingRadius > 0)
							{
								for (i = 0; i < this->LiftingRadius; i++)
								{
									this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
									TempList = VerticalList;
									VerticalList = HorizontalList;
									HorizontalList = TempList;
								}
							}

							// all v in k-disk are added to R^j
							int Range = HorizontalList->GetNumberOfIds();
							for (i = 0; i < Range; i++)
							{
								c = HorizontalList->GetId(i);
								if (c < NOldPoints)
									R_next[HorizontalList->GetId(i)] = true;
							}
						}
					}

					if (this->ArithmeticType == 1)
					{
						// recalculate w in hat(W)^j
						vtkIntArray *IntWavelets = this->Filters[j]->IntegerWavelets;
						int *IntWavelet;
						double P1[3], P2[3], V1[3], V[3];
						for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
						{
							if (W_next[PtId - NOldPoints] == true)
							{
								p1 = this->Filters[j]->vertices[PtId].parent1;
								p2 = this->Filters[j]->vertices[PtId].parent2;
								NewPoints->GetPoint(PtId, V1);
								NewPoints->GetPoint(p1, P1);
								NewPoints->GetPoint(p2, P2);

								IntWavelet = IntWavelets->GetPointer((PtId - NOldPoints) * 3);

								IntWavelet[0] = (int)floor(V1[0] - 0.5*(P1[0] + P2[0]));
								IntWavelet[1] = (int)floor(V1[1] - 0.5*(P1[1] + P2[1]));
								IntWavelet[2] = (int)floor(V1[2] - 0.5*(P1[2] + P2[2]));
							}
						}
						IntWavelets->Modified();

						// recalulate v in R^j
						for (PtId = 0; PtId < NOldPoints; PtId++)
						{
							if (R_next[PtId] == true)
							{
								V[0] = 0;
								V[1] = 0;
								V[2] = 0;

								NonZero = this->Filters[j]->Alpha->FirstHor[PtId];

								while (NonZero)
								{
									IntWavelet = IntWavelets->GetPointer((NonZero->j) * 3);

									V[0] += NonZero->Value*IntWavelet[0];
									V[1] += NonZero->Value*IntWavelet[1];
									V[2] += NonZero->Value*IntWavelet[2];

									NonZero = NonZero->NextHor;
								}

								NewPoints->GetPoint(PtId, P1);

								P2[0] = P1[0] + floor(V[0] + 0.5);
								P2[1] = P1[1] + floor(V[1] + 0.5);
								P2[2] = P1[2] + floor(V[2] + 0.5);

								OldPoints->SetPoint(PtId, P2);
							}
						}
					}

					R_raw = R_next;
				}
				end = clock();

				out.open(outinfo, ios::app);
				out << double(end - start) / CLOCKS_PER_SEC;
				out.close();
			}

			///////////////////////////////////////////////////////////
			////////// count modified vertices and wavelets ///////////
			///////////////////////////////////////////////////////////
			{
				// number of vertices and wavelets to be modified
				int number_of_modified_wavelets = 0;
				int number_of_modified_vertices = 0;

				// R^finest initialize
				bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
				for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
				{
					R_raw[PtId] = ModifiedTable[PtId];
					if (R_raw[PtId] == true)
						number_of_modified_vertices += 1;
				}

				vtkIdList *Edges = vtkIdList::New();
				vtkIdType PtId;
				vtkIdType edge;
				vtkIdType c, p1, p2;

				vtkIdList *HorizontalList, *VerticalList, *TempList;
				HorizontalList = vtkIdList::New();
				VerticalList = vtkIdList::New();

				vtkSparseMatrix::NonZeroElement *NonZero;

				for (j = 0; j < this->NumberOfFilters; j++)
				{
					vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
					vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
					vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
					vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

					bool *R_next = new bool[NOldPoints]();
					bool *W_next = new bool[NNewPoints - NOldPoints]();

					// w impacted by v in R^(j+1) must be added to hat(W)^j
					//     and propagated into k-disk of its parents
					for (PtId = 0; PtId < NOldPoints; PtId++)
					{
						if (R_raw[PtId] == true)
						{
							// add v to R^j 
							R_next[PtId] = true;

							this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
							for (e = 0; e < Edges->GetNumberOfIds(); e++)
							{
								// find c impacted by v
								edge = Edges->GetId(e);
								c = this->Filters[j]->edgesvector[edge].child;
								if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
								{
									// add c to hat(W)^j
									W_next[c - NOldPoints] = true;

									// propagate parents of c into k-disk
									HorizontalList->SetNumberOfIds(2);
									HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
									HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

									if (this->LiftingRadius > 0)
									{
										for (i = 0; i < this->LiftingRadius; i++)
										{
											this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
											TempList = VerticalList;
											VerticalList = HorizontalList;
											HorizontalList = TempList;
										}
									}

									// all v in k-disk are added to R^j
									int Range = HorizontalList->GetNumberOfIds();
									for (i = 0; i < Range; i++)
									{
										c = HorizontalList->GetId(i);
										if (c < NOldPoints)
										{
											R_next[HorizontalList->GetId(i)] = true;
										}
									}
								}
							}
						}
					}
					// mergable v in R^(j+1) must be added to hat(W)^j
					for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
					{
						if (R_raw[PtId] == true)
						{
							W_next[PtId - NOldPoints] = true;

							// propagate parents of PtId into k-disk
							HorizontalList->SetNumberOfIds(2);
							HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
							HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

							if (this->LiftingRadius > 0)
							{
								for (i = 0; i < this->LiftingRadius; i++)
								{
									this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
									TempList = VerticalList;
									VerticalList = HorizontalList;
									HorizontalList = TempList;
								}
							}

							// all v in k-disk are added to R^j
							int Range = HorizontalList->GetNumberOfIds();
							for (i = 0; i < Range; i++)
							{
								c = HorizontalList->GetId(i);
								if (c < NOldPoints)
								{
									R_next[HorizontalList->GetId(i)] = true;
								}
							}
						}
					}

					for (PtId = 0; PtId < NOldPoints; PtId++)
					{
						if (R_next[PtId] == true)
							number_of_modified_vertices += 1;
					}
					for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
					{
						if (W_next[PtId - NOldPoints] == true)
							number_of_modified_wavelets += 1;
					}

					R_raw = R_next;
				}

				out.open(outinfo, ios::app);
				out << " " << t;
				out << " " << number_of_modified_vertices + number_of_modified_wavelets << endl;
				out.close();
			}
		}
	}
}

// modify k% of entire vertices by propagating to neighbors (if not [k,k+0.5]%, rechoose random initial point )
void vtkMultiresolutionIO::RoiLocal(int k)
{
	int i, j, t, e;
	vtkIdType PtId;

	std::string outinfo = this->InputFileName;
	outinfo += "_info.txt";
	std::ofstream out;

	out.open(outinfo, ios::app);
	out << "ROI: ";
	out.close();

	this->MeshWindow->SetInputData(this->SynthesisMeshes[0]);

	int number_of_vertices = this->Filters[0]->GetOutput()->GetNumberOfPoints();
	int number_of_modified = 0;

	// random setting
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, number_of_vertices-1); // define the range
	int PtIdR;

	// choose modified local ROI vertices (init with all 0s)
	bool *ModifiedTable = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
	while (1)
	{
		// choose random initial point
		PtIdR = distr(eng);

		vtkIdList *HorizontalList, *VerticalList, *TempList;
		HorizontalList = vtkIdList::New();
		VerticalList = vtkIdList::New();

		// propagate PtIdR into about k% of entire vertices
		HorizontalList->SetNumberOfIds(1);
		HorizontalList->SetId(0, PtIdR);

		if (this->LiftingRadius >= 0)
		{
			bool isolated = false;
			while ( HorizontalList->GetNumberOfIds() < (0.01 * k * number_of_vertices) )
			{
				this->Filters[0]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
				if (HorizontalList->GetNumberOfIds() == VerticalList->GetNumberOfIds())
				{
					isolated = true;
					break;
				}
				TempList = VerticalList;
				VerticalList = HorizontalList;
				HorizontalList = TempList;
			}

			if (isolated == true)
				continue;
			
			// check ROI size < (k+0.5)%
			if (HorizontalList->GetNumberOfIds() < (0.01 * (k + 0.5) * number_of_vertices))
			{
				// all v in local ROI are added to ModifiedTable
				number_of_modified = HorizontalList->GetNumberOfIds();
				for (i = 0; i < number_of_modified; i++)
					ModifiedTable[HorizontalList->GetId(i)] = true;

				break;
			}
		}

		cout << "Failed " << PtIdR << " vertex index" << endl;
	}

	clock_t start, end;

	// 2) with lifting
	//      after approximate,
	//		v is recalculated from Mesh2
	//		all vertices must be recalculated sequentially from finest to coarsest
	//      all SynthesisMeshes must be handled.
	if (this->Lifting == 1)
	{
		// propagation //
		{
			start = clock();

			// R^finest initialize
			bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
			for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
				R_raw[PtId] = ModifiedTable[PtId];

			//// R^finest modify ///////////////////////////////////// move 0.1
			//double PF[3];
			//double Bounds[6];
			//double xdyn, ydyn, zdyn;
			//this->Filters[0]->GetOutput()->GetPoints()->GetBounds(Bounds);
			//xdyn = Bounds[1] - Bounds[0];
			//ydyn = Bounds[3] - Bounds[2];
			//zdyn = Bounds[5] - Bounds[4];
			//for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
			//{
			//	if (R_raw[PtId] == true)
			//	{
			//		this->Filters[0]->GetOutput()->GetPoints()->GetPoint(PtId, PF);

			//		PF[0] = PF[0] + xdyn*0.2;
			//		PF[1] = PF[1] + ydyn*0.2;
			//		PF[2] = PF[2] + zdyn*0.2;

			//		this->Filters[0]->GetOutput()->GetPoints()->SetPoint(PtId, PF);
			//	}
			//}

			vtkIdList *Edges = vtkIdList::New();
			vtkIdType PtId;
			vtkIdType edge;
			vtkIdType c, p1, p2;

			vtkIdList *HorizontalList, *VerticalList, *TempList;
			HorizontalList = vtkIdList::New();
			VerticalList = vtkIdList::New();

			vtkSparseMatrix::NonZeroElement *NonZero;

			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
				vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
				vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
				vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

				bool *R_next = new bool[NOldPoints]();
				bool *W_next = new bool[NNewPoints - NOldPoints]();

				int lengh_of_R_curr = 0;
				
				//////////////////////////////////////
				////  4 step rules - (a) and (b)  ////
				//////////////////////////////////////
				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						// add v to R^j 
						R_next[PtId] = true;

						this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
						for (e = 0; e < Edges->GetNumberOfIds(); e++)
						{
							// find c impacted by v
							edge = Edges->GetId(e);
							c = this->Filters[j]->edgesvector[edge].child;
							if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
							{
								// add c to hat(W)^j
								W_next[c - NOldPoints] = true; //////////////////////////////////////// cannot enter this???

								// propagate parents of c into k-disk
								HorizontalList->SetNumberOfIds(2);
								HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
								HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

								if (this->LiftingRadius > 0)
								{
									for (i = 0; i < this->LiftingRadius; i++)
									{
										this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
										TempList = VerticalList;
										VerticalList = HorizontalList;
										HorizontalList = TempList;
									}
								}

								// all v in k-disk are added to R^j
								int Range = HorizontalList->GetNumberOfIds();
								for (i = 0; i < Range; i++)
								{
									c = HorizontalList->GetId(i);
									if (c < NOldPoints)
										R_next[HorizontalList->GetId(i)] = true;
								}
							}
						}
					}
				}
				//////////////////////////////////////
				////  4 step rules - (c) and (d)  ////
				//////////////////////////////////////
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						// mergable v in R^(j+1) must be added to hat(W)^j
						W_next[PtId - NOldPoints] = true;

						// propagate parents of PtId into k-disk
						HorizontalList->SetNumberOfIds(2);
						HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
						HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

						if (this->LiftingRadius > 0)
						{
							for (i = 0; i < this->LiftingRadius; i++)
							{
								this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
								TempList = VerticalList;
								VerticalList = HorizontalList;
								HorizontalList = TempList;
							}
						}

						// all v in k-disk are added to R^j
						int Range = HorizontalList->GetNumberOfIds();
						for (i = 0; i < Range; i++)
						{
							c = HorizontalList->GetId(i);
							if (c < NOldPoints)
								R_next[HorizontalList->GetId(i)] = true;
						}
					}
				}

				//// display before recalculate
				//this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
				//if (this->Display == 2)
				//{
				//	// test. set face colors
				//	// fail... how...???
				//	vtkIntArray* colors = vtkIntArray::New();
				//	colors->SetNumberOfValues(NNewPoints);
				//	for (i = 0; i < NNewPoints; i++)
				//	{
				//		if (R_raw[i] == true)
				//		{
				//			colors->SetValue(i, 1);
				//		}
				//		else
				//		{
				//			colors->SetValue(i, 0);
				//		}
				//	}
				//	double c1[3] = { 0.8,0.9,0 }; // yellow
				//	double c2[3] = { 255. / 255, 255. / 255, 255. / 255 }; // white
				//	this->MeshWindow->DisplayVerticesColors2(colors, c1, c2);
				//	this->MeshWindow->Render();
				//	this->MeshWindow->Interact();
				//}

				if (this->ArithmeticType == 1)
				{
					// recalculate w in hat(W)^j
					vtkIntArray *IntWavelets = this->Filters[j]->IntegerWavelets;
					int *IntWavelet;
					double P1[3], P2[3], V1[3], V[3];
					for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
					{
						if (W_next[PtId - NOldPoints] == true)
						{
							p1 = this->Filters[j]->vertices[PtId].parent1;
							p2 = this->Filters[j]->vertices[PtId].parent2;
							NewPoints->GetPoint(PtId, V1);
							NewPoints->GetPoint(p1, P1);
							NewPoints->GetPoint(p2, P2);

							IntWavelet = IntWavelets->GetPointer((PtId - NOldPoints) * 3);

							IntWavelet[0] = (int)floor(V1[0] - 0.5*(P1[0] + P2[0]));
							IntWavelet[1] = (int)floor(V1[1] - 0.5*(P1[1] + P2[1]));
							IntWavelet[2] = (int)floor(V1[2] - 0.5*(P1[2] + P2[2]));
						}
					}
					IntWavelets->Modified();

					// recalulate v in R^j
					for (PtId = 0; PtId < NOldPoints; PtId++)
					{
						if (R_next[PtId] == true)
						{
							V[0] = 0;
							V[1] = 0;
							V[2] = 0;

							NonZero = this->Filters[j]->Alpha->FirstHor[PtId];

							while (NonZero)
							{
								IntWavelet = IntWavelets->GetPointer((NonZero->j) * 3);

								V[0] += NonZero->Value*IntWavelet[0];
								V[1] += NonZero->Value*IntWavelet[1];
								V[2] += NonZero->Value*IntWavelet[2];

								NonZero = NonZero->NextHor;
							}

							NewPoints->GetPoint(PtId, P1);

							P2[0] = P1[0] + floor(V[0] + 0.5);
							P2[1] = P1[1] + floor(V[1] + 0.5);
							P2[2] = P1[2] + floor(V[2] + 0.5);

							OldPoints->SetPoint(PtId, P2);
						}
					}
				}

				this->MeshWindow->SetInputData(this->Filters[j]->GetSubdivisionInput());
				if (this->Display == 2)
				{
					// test. set face colors
					// fail... how...???
					vtkIntArray* colors = vtkIntArray::New();
					colors->SetNumberOfValues(NOldPoints);
					for (i = 0; i < NOldPoints; i++)
					{
						if (R_next[i] == true)
						{
							colors->SetValue(i, 1);
						}
						else
						{
							colors->SetValue(i, 0);
						}
					}
					double c1[3] = { 0.8,0.9,0 }; // yellow
					double c2[3] = { 60. / 255, 63. / 255, 200. / 255 }; // blue
					this->MeshWindow->DisplayVerticesColors2(colors, c1, c2);
					this->MeshWindow->Render();
					this->MeshWindow->Interact();
				}

				this->MeshWindow->SetInputData(this->Filters[j]->GetOutput());
				if (this->Display == 2)
				{
					vtkIntArray* colors = vtkIntArray::New();
					colors->SetNumberOfValues(NNewPoints);
					for (i = 0; i < NOldPoints; i++)
					{
						colors->SetValue(i, 0);
					}
					for (i = NOldPoints; i < NNewPoints; i++)
					{
						if (W_next[i - NOldPoints] == true)
						{
							colors->SetValue(i, 1);
						}
						else
						{
							colors->SetValue(i, 0);
						}
					}
					double c1[3] = { 0.8,0.9,0 }; // yellow
					double c2[3] = { 0. / 255, 244. / 255, 1. / 255 }; // green
					this->MeshWindow->DisplayVerticesColors2(colors, c1, c2);
					this->MeshWindow->Render();
					this->MeshWindow->Interact();
				}

				R_raw = R_next;
			}
			end = clock();

			out.open(outinfo, ios::app);
			out << double(end - start) / CLOCKS_PER_SEC << endl;
			out.close();
		}

		///////////////////////////////////////////////////////////
		////////// count modified vertices and wavelets ///////////
		///////////////////////////////////////////////////////////
		{
			// number of vertices and wavelets to be modified
			int number_of_modified_wavelets = 0;
			int number_of_modified_vertices = 0;

			// R^finest initialize
			bool *R_raw = new bool[this->Filters[0]->GetOutput()->GetNumberOfPoints()]();
			for (PtId = 0; PtId < this->Filters[0]->GetOutput()->GetNumberOfPoints(); PtId++)
			{
				R_raw[PtId] = ModifiedTable[PtId];
				if (R_raw[PtId] == true)
					number_of_modified_vertices += 1;
			}

			vtkIdList *Edges = vtkIdList::New();
			vtkIdType PtId;
			vtkIdType edge;
			vtkIdType c, p1, p2;

			vtkIdList *HorizontalList, *VerticalList, *TempList;
			HorizontalList = vtkIdList::New();
			VerticalList = vtkIdList::New();

			vtkSparseMatrix::NonZeroElement *NonZero;

			for (j = 0; j < this->NumberOfFilters; j++)
			{
				vtkIdType NOldPoints = this->Filters[j]->GetSubdivisionInput()->GetNumberOfPoints();
				vtkIdType NNewPoints = this->Filters[j]->GetOutput()->GetNumberOfPoints();
				vtkPoints *OldPoints = this->Filters[j]->GetSubdivisionInput()->GetPoints();
				vtkPoints *NewPoints = this->Filters[j]->GetOutput()->GetPoints();

				bool *R_next = new bool[NOldPoints]();
				bool *W_next = new bool[NNewPoints - NOldPoints]();

				// w impacted by v in R^(j+1) must be added to hat(W)^j
				//     and propagated into k-disk of its parents
				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						// add v to R^j 
						R_next[PtId] = true;

						this->Filters[j]->GetSubdivisionInput()->GetVertexNeighbourEdges(PtId, Edges);
						for (e = 0; e < Edges->GetNumberOfIds(); e++)
						{
							// find c impacted by v
							edge = Edges->GetId(e);
							c = this->Filters[j]->edgesvector[edge].child;
							if ((c >= NOldPoints) && (this->Filters[j]->GetOutput()->IsEdge(PtId, c) >= 0))
							{
								// add c to hat(W)^j
								W_next[c - NOldPoints] = true;

								// propagate parents of c into k-disk
								HorizontalList->SetNumberOfIds(2);
								HorizontalList->SetId(0, this->Filters[j]->vertices[c].parent1);
								HorizontalList->SetId(1, this->Filters[j]->vertices[c].parent2);

								if (this->LiftingRadius > 0)
								{
									for (i = 0; i < this->LiftingRadius; i++)
									{
										this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
										TempList = VerticalList;
										VerticalList = HorizontalList;
										HorizontalList = TempList;
									}
								}

								// all v in k-disk are added to R^j
								int Range = HorizontalList->GetNumberOfIds();
								for (i = 0; i < Range; i++)
								{
									c = HorizontalList->GetId(i);
									if (c < NOldPoints)
									{
										R_next[HorizontalList->GetId(i)] = true;
									}
								}
							}
						}
					}
				}
				// mergable v in R^(j+1) must be added to hat(W)^j
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (R_raw[PtId] == true)
					{
						W_next[PtId - NOldPoints] = true;

						// propagate parents of PtId into k-disk
						HorizontalList->SetNumberOfIds(2);
						HorizontalList->SetId(0, this->Filters[j]->vertices[PtId].parent1);
						HorizontalList->SetId(1, this->Filters[j]->vertices[PtId].parent2);

						if (this->LiftingRadius > 0)
						{
							for (i = 0; i < this->LiftingRadius; i++)
							{
								this->Filters[j]->GetSubdivisionInput()->GetNeighbours(HorizontalList, VerticalList);
								TempList = VerticalList;
								VerticalList = HorizontalList;
								HorizontalList = TempList;
							}
						}

						// all v in k-disk are added to R^j
						int Range = HorizontalList->GetNumberOfIds();
						for (i = 0; i < Range; i++)
						{
							c = HorizontalList->GetId(i);
							if (c < NOldPoints)
							{
								R_next[HorizontalList->GetId(i)] = true;
							}
						}
					}
				}

				for (PtId = 0; PtId < NOldPoints; PtId++)
				{
					if (R_next[PtId] == true)
						number_of_modified_vertices += 1;
				}
				for (PtId = NOldPoints; PtId < NNewPoints; PtId++)
				{
					if (W_next[PtId - NOldPoints] == true)
						number_of_modified_wavelets += 1;
				}

				R_raw = R_next;
			}

			out.open(outinfo, ios::app);
			out << "started at " << PtIdR << endl;
			out << " " << number_of_modified;
			out << " / " << number_of_vertices;
			out << " = " << (double)number_of_modified / number_of_vertices << "%" << endl;
			out << " " << number_of_modified_vertices;
			out << " " << number_of_modified_wavelets;
			out << " " << number_of_modified_vertices + number_of_modified_wavelets << endl;
			out.close();
		}
	}
}

void vtkMultiresolutionIO::SetInput2()
{
	this->MeshWindow->SetInputData(this->Input);
}

void vtkMultiresolutionIO::ExportOBJ(vtkSurface *Surface)
{
	// write finest mesh to file
	if (this->WriteOutput > 0)
	{
		std::stringstream strfile;

		if (this->FileType == 0)
		{
			std::string str1 = this->InputFileName;
			str1 += "_mesh.iv";

			strfile << str1;
			Surface->WriteInventor(strfile.str().c_str());
		}
		else
		{
			std::string str1 = this->InputFileName;
			str1 += "_mesh.obj";

			strfile << str1;
			std::ofstream out(strfile.str());

			int NumOfPoints = Surface->GetNumberOfPoints();
			vtkPoints *Points = Surface->GetPoints();
			double p[3];

			if (this->ArithmeticType == 1)
			{
				double Factor, Tx, Ty, Tz;
				this->Filters[0]->GetOutput()->GetScalingFactors(Factor, Tx, Ty, Tz);

				for (int i = 0; i < NumOfPoints; i++)
				{
					Points->GetPoint(i, p);
					p[0] = (p[0] / Factor) + Tx;
					p[1] = (p[1] / Factor) + Ty;
					p[2] = (p[2] / Factor) + Tz;
					out << "v" << " " << p[0] << " " << p[1] << " " << p[2] << endl;
				}
			}
			else
			{
				for (int i = 0; i < NumOfPoints; i++)
				{
					Points->GetPoint(i, p);
					out << "v" << " " << p[0] << " " << p[1] << " " << p[2] << endl;
				}
			}

			int NumOfFaces = Surface->GetNumberOfCells();

			for (int i = 0; i < NumOfFaces; i++)
			{
				vtkIdType v1, v2, v3;
				Surface->GetFaceVertices(i,v1,v2,v3);
				out << "f" << " " << v1+1 << " " << v2+1 << " " << v3+1 << endl;
			}

			out.close();
		}

		if (this->ArithmeticType == 1)
			Surface->QuantizeCoordinates(this->Quantization);
	}
}


vtkMultiresolutionIO::vtkMultiresolutionIO()
{
	this->EdgeAngleThreshold=0.3;
	this->GeometricConstraint=0;
	this->NumberOfFilters=0;
	this->WGC=0.25;
	this->Quantization=12;
	this->ArithmeticType=1;
	this->Lifting=0;
	this->LiftingRadius=1;
	this->DisplayEfficiency=0;
	this->Display=0;
	this->DisplayText=0;
	this->WriteRepport=1;
	this->NumberOfBitPlanes=12;
	this->NumberOfStartBitPlanes=1;
	this->WriteOutput=1;  // -1: print all meshes, n>0: print n finest mesh only
	this->GeometryPrediction=0;
	this->MaxNumberOfLevels=99;
	this->Input=0;
	this->Output=0;
	this->FileType=1;

	this->MeshWindow=RenderWindow::New();
	this->PointsIds=vtkIdList::New();
	this->SignificanceCodes=0;
	this->Timer=vtkTimerLog::New();
	this->InputFileName = NULL;
	this->SetFileName("out.ddd");
}

vtkMultiresolutionIO::~vtkMultiresolutionIO()
{
	this->PointsIds->Delete();
	int i;
	for (i=0;i<this->NumberOfFilters;i++)
		this->Filters[i]->Delete();
	
	if (this->NumberOfFilters!=0)
		this->SynthesisMeshes[this->NumberOfFilters]->Delete();
		
	if (this->Input)
		this->Input->UnRegister(this);

	if (this->Output)
		this->Output->UnRegister(this);

	this->MeshWindow->Delete();
	this->Timer->Delete();
}
