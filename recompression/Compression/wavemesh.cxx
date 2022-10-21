/*=========================================================================

  Program:   Wavemesh : a progressive lossless compression scheme for 3D triagular meshes
  Module:    wavemesh.cxx
  Language:  C++
  Date:      2008/08
  Auteur:    Sebastien Valette
    This software is governed by the GPL license (see License.txt)
=========================================================================*/
// .NAME wavemesh 
// .SECTION Description
/*==========================================================================

  Date:      2021/10
  Author:    Yeonghun Kim
  
============================================================================*/

#include "vtkMultiresolutionIO.h"
#include <time.h>
#include <string>
#include <iostream>
#include <filesystem>

int main( int argc, char *argv[] )
{

	if(argc<2)
	{
		cout<<"Usage :"<<endl;
		cout<<"Compression : wavemesh c file [options]"<<endl;
		cout<<"Decompression : wavemesh d file [options]"<<endl;
		cout << "Test for SWStarLab 2022 : wavemesh test file [options]" << endl;
		cout<<endl;
		cout<<"Optionnal arguments :"<<endl;
		cout<<"-d 0/1/2 changes display type (default : 0):"<<endl;
		cout<<"              0: no display"<<endl;
		cout<<"              1: display all levels"<<endl;
		cout<<"              2: display with stop at each level"<<endl;
		cout<<"-a 0/1 changes arithmetics : 0=floats (zerotree coding)"<<endl;
		cout<<"                             1=integers (original wavemesh coding)"<<endl;
		cout<<"-q n : defines the coordinates quantization to n bits (default : 12 bits)"<<endl;
		cout<<"-b n : defines the number of bitplanes for zerotree coding (default : 12)"<<endl;
		cout<<"-l r : defines the lifting radius (-1=no lifting. Default : -1)"<<endl;
		cout<<"-g 0/1 : enables/disable geometric constraints for simplification"<<endl;
		cout<<"-et threshold : defines the edge angle threshold (default : 0.3)"<<endl;
		cout<<"-wt threshold : defines the wavelet ratio threshold (default : 0.25)"<<endl;
		cout<<"-o filename : defines the compressed output file name (default : .out)"<<endl;
		
		return (0);
	}
	else
	{
		//cout << "For mor help on wavemesh, execute the program without any argument" << endl;
	}

	namespace fs = std::filesystem;

	std::string path = argv[2];
	vtkMultiresolutionIO *MIO = vtkMultiresolutionIO::New();

	// input file name setting
	std::string inputfile = path;
	inputfile = inputfile.substr(0, inputfile.size() - 4);
	MIO->SetInputFileName(inputfile.c_str());

	int Arithmetics = 1;

	MIO->SetLifting(1);
	MIO->SetLiftingRadius(2);
	MIO->SetQuantization(24);

	// Parse optionnal arguments
	int ArgumentsIndex = 3;
	while (ArgumentsIndex < argc)
	{
		if (strcmp(argv[ArgumentsIndex], "-d") == 0)
		{
			int Display = atoi(argv[ArgumentsIndex + 1]);
			MIO->SetDisplay(Display);
			cout << "Display=" << Display << endl;
		}

		if (strcmp(argv[ArgumentsIndex], "-a") == 0)
		{
			Arithmetics = atoi(argv[ArgumentsIndex + 1]);
			cout << "Arithmetics (0=floats (zerotree) ; 1=integers (original wavemesh)): " << Arithmetics << endl;
			MIO->SetArithmeticType(Arithmetics);
		}

		if (strcmp(argv[ArgumentsIndex], "-q") == 0)
		{
			int Quantization = atoi(argv[ArgumentsIndex + 1]);
			//cout << "Coordinates Quantization :" << Quantization << endl;
			MIO->SetQuantization(Quantization);
		}

		if (strcmp(argv[ArgumentsIndex], "-b") == 0)
		{
			int NumberOfBitPlanes = atoi(argv[ArgumentsIndex + 1]);
			cout << "Number Of BitPlanes :" << NumberOfBitPlanes << endl;
			MIO->SetNumberOfBitPlanes(NumberOfBitPlanes);
		}

		if (strcmp(argv[ArgumentsIndex], "-l") == 0)
		{
			int Lifting = atoi(argv[ArgumentsIndex + 1]);
			if (Lifting == -1)
			{
				cout << "Lifting scheme deactivated" << endl;
				MIO->SetLifting(0);
			}
			else
			{
				MIO->SetLifting(1);
				MIO->SetLiftingRadius(Lifting);
				// cout << "Lifting radius : " << Lifting << endl;
			}
		}

		if (strcmp(argv[ArgumentsIndex], "-g") == 0)
		{
			int Geometry = atoi(argv[ArgumentsIndex + 1]);
			cout << "Geometric constraints for simplification :" << Geometry << endl;
			MIO->SetGeometricalConstraint(Geometry);
		}

		if (strcmp(argv[ArgumentsIndex], "-et") == 0)
		{
			double Threshold = atof(argv[ArgumentsIndex + 1]);
			cout << "Edge angle threshold :" << Threshold << endl;
			MIO->SetEdgeAngleThreshold(Threshold);
		}

		if (strcmp(argv[ArgumentsIndex], "-wt") == 0)
		{
			double Threshold = atof(argv[ArgumentsIndex + 1]);
			cout << "Wavelet ratio threshold :" << Threshold << endl;
			MIO->SetWGC(Threshold);
		}

		if (strcmp(argv[ArgumentsIndex], "-o") == 0)
		{
			cout << "Output File : " << argv[ArgumentsIndex + 1] << endl;
			MIO->SetFileName(argv[ArgumentsIndex + 1]);
		}
		ArgumentsIndex += 2;
	}

	if (strcmp(argv[1], "c") == 0)
	{
		vtkSurface *Mesh = vtkSurface::New();
		cout << "Load : " << path << endl;
		Mesh->CreateFromFile(path.c_str());
		Mesh->DisplayMeshProperties();

		if (Arithmetics == 1)
			Mesh->QuantizeCoordinates(MIO->GetQuantization());

		MIO->SetInput(Mesh);

		// output file name (instead of '.out')
		std::string strout = MIO->GetInputFileName();
		strout += ".out";
		MIO->SetFileName(strout.c_str());

		MIO->Analyse();
		MIO->Synthetize();
		MIO->Approximate();
		MIO->Write();
		Mesh->Delete();
		MIO->Delete();
	}

	// Compress -> Modify it (from _m.obj file) -> Recompress -> Decompress
	if (strcmp(argv[1], "r") == 0)
	{
		// merge first (since it has no parent-child relations)
		vtkSurface *Mesh = vtkSurface::New();
		cout << "Load : " << path << endl;
		Mesh->CreateFromFile(path.c_str());
		Mesh->DisplayMeshProperties();

		std::string strout = MIO->GetInputFileName();
		MIO->SetFileName(strout.c_str());

		if (Arithmetics == 1)
			Mesh->QuantizeCoordinates(MIO->GetQuantization());

		MIO->SetInput(Mesh);

		MIO->Analyse();  // you must not move vertices in this step!
		MIO->Synthetize();

		clock_t start;
		clock_t end;

		if (MIO->GetLifting())
		{
			// re-compression
			while (1)
			{
				MIO->RoiSetting();

				// timing experiment
				start = clock();
				MIO->Approximate();
				end = clock();
				cout << "Approximate elapsed time: " << (double)(end - start) << endl;

				// timing experiment
				start = time(NULL);
				MIO->Write();
				end = time(NULL);
				cout << "Write elapsed time: " << (double)(end - start) << endl;
				break;
			}
		}
		else
		{
			MIO->Approximate();
			// re-compression
			while (1)
			{
				MIO->RoiSetting();

				// timing experiment
				start = clock();
				MIO->Write();
				end = clock();
				cout << "Write elapsed time: " << (double)(end - start) << endl;
				break;
			}
		}

		Mesh->Delete();
		MIO->Delete();

		// test 
		MIO = vtkMultiresolutionIO::New();

		// Parse optionnal arguments
		ArgumentsIndex = 3;
		while (ArgumentsIndex < argc)
		{
			if (strcmp(argv[ArgumentsIndex], "-d") == 0)
			{
				int Display = atoi(argv[ArgumentsIndex + 1]);
				MIO->SetDisplay(Display);
				cout << "Display=" << Display << endl;
			}

			if (strcmp(argv[ArgumentsIndex], "-a") == 0)
			{
				Arithmetics = atoi(argv[ArgumentsIndex + 1]);
				cout << "Arithmetics (0=floats (zerotree) ; 1=integers (original wavemesh)): " << Arithmetics << endl;
				MIO->SetArithmeticType(Arithmetics);
			}

			if (strcmp(argv[ArgumentsIndex], "-q") == 0)
			{
				int Quantization = atoi(argv[ArgumentsIndex + 1]);
				cout << "Coordinates Quantization :" << Quantization << endl;
				MIO->SetQuantization(Quantization);
			}

			if (strcmp(argv[ArgumentsIndex], "-b") == 0)
			{
				int NumberOfBitPlanes = atoi(argv[ArgumentsIndex + 1]);
				cout << "Number Of BitPlanes :" << NumberOfBitPlanes << endl;
				MIO->SetNumberOfBitPlanes(NumberOfBitPlanes);
			}

			if (strcmp(argv[ArgumentsIndex], "-l") == 0)
			{
				int Lifting = atoi(argv[ArgumentsIndex + 1]);
				if (Lifting == -1)
				{
					cout << "Lifting scheme deactivated" << endl;
					MIO->SetLifting(0);
				}
				else
				{
					MIO->SetLifting(1);
					MIO->SetLiftingRadius(Lifting);
					// cout << "Lifting radius : " << Lifting << endl;
				}
			}

			if (strcmp(argv[ArgumentsIndex], "-g") == 0)
			{
				int Geometry = atoi(argv[ArgumentsIndex + 1]);
				cout << "Geometric constraints for simplification :" << Geometry << endl;
				MIO->SetGeometricalConstraint(Geometry);
			}

			if (strcmp(argv[ArgumentsIndex], "-et") == 0)
			{
				double Threshold = atof(argv[ArgumentsIndex + 1]);
				cout << "Edge angle threshold :" << Threshold << endl;
				MIO->SetEdgeAngleThreshold(Threshold);
			}

			if (strcmp(argv[ArgumentsIndex], "-wt") == 0)
			{
				double Threshold = atof(argv[ArgumentsIndex + 1]);
				cout << "Wavelet ratio threshold :" << Threshold << endl;
				MIO->SetWGC(Threshold);
			}

			if (strcmp(argv[ArgumentsIndex], "-o") == 0)
			{
				cout << "Output File : " << argv[ArgumentsIndex + 1] << endl;
				MIO->SetFileName(argv[ArgumentsIndex + 1]);
			}
			ArgumentsIndex += 2;
		}

		// check
		MIO->Read();  // .out
		MIO->Delete();
	}

	// d : decompression -> mesh.ply ( WriteOutput = -1 (all mesh) or n>0 (n finest mesh) )
	if (strcmp(argv[1], "d") == 0)
	{
		std::string strout = path;
		MIO->SetFileName(strout.c_str());

		clock_t start, end;
		start = clock();
		MIO->Read();
		end = clock();

		int write_info = 1;

		if (write_info == 1)
		{
			std::string infoout = MIO->GetInputFileName();
			infoout += "_decoding.txt";

			std::ofstream out;
			out.open(infoout, ios::trunc);
			out << "Decoding elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();
		}

		MIO->Delete();
	}


	// test :
	// 1. compression -> mesh.out (print Analyse, Synthetize, Approximate, Write, Total Compression Time)
	// 2. decompression -> mesh_d.ply (print Total Decompression Time (= Read))
	// 3. modification -> mesh_m.ply (print 'randomly modify all geometry')
	// 4. recompression -> mesh_m.out (print Total Recompression Time)
	// //5. decompression -> mesh_md.out (need...?)
	if (strcmp(argv[1], "test") == 0)
	{
		// Mesh loading //
			
		std::string basename = MIO->GetInputFileName();
		std::string strout = MIO->GetInputFileName(); strout += ".out";

		clock_t start, end;
		double total_elapsed = 0.;
		double re_elapsed = 0.;

		vtkSurface* Mesh = vtkSurface::New();
		Mesh->CreateFromFile(path.c_str());

		cout << "Loaded : " << path << endl;
		cout << "File size: " << (double)fs::file_size(path) / (1024 * 1024) << " MB" << endl;
		cout << endl;

		// Compression //

		cout << "Compression in progress ... " << endl;

		if (Arithmetics == 1)
			Mesh->QuantizeCoordinates(MIO->GetQuantization());

		MIO->SetInput(Mesh);
		MIO->SetFileName(strout.c_str());	// output file name (instead of '.out')

		//start = clock();
		MIO->Analyse();
		//cout << "Analyse. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		//start = clock();
		MIO->Synthetize();
		//cout << "Synthetize. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		//start = clock();
		MIO->Approximate();
		//cout << "Approximate. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		//start = clock();
		MIO->Write();
		//cout << "Write. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		// Mesh->Delete();

		cout << "completed" << endl;

		// Decompression //

		cout << "Decompression in progress ... "  << endl;

		vtkMultiresolutionIO* MIO_d = vtkMultiresolutionIO::New();

		start = clock();
		MIO_d->SetInputFileName(basename.c_str()); // input file name setting
		MIO_d->SetFileName(strout.c_str());
		MIO_d->Read();
		end = clock();
		MIO_d->Delete();

		cout << "completed. " << (double)(end - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		// Recompression //

		cout << "Recompression in progress ... "  << endl;

		start = clock();
		std::string strnewout = MIO->GetInputFileName(); strnewout += "_m.out";
		MIO->SetFileName(strnewout.c_str());
		MIO->RoiRandom(1);
		MIO->Write();
		end = clock();
		MIO->Delete();
		Mesh->Delete();

		cout << "completed. " << (double)(end - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;
	}

	
	// test for 2022 test :
	if (strcmp(argv[1], "testdev") == 0)
	{
		// Mesh loading //

		std::string basename = MIO->GetInputFileName();
		std::string strout = MIO->GetInputFileName(); strout += ".out";


		///////////////// DO NOT NEED IF COMPRESSION END
		clock_t start, end;
		double total_elapsed = 0.;
		double re_elapsed = 0.;

		vtkSurface* Mesh = vtkSurface::New();
		Mesh->CreateFromFile(path.c_str());

		cout << "Loaded : " << path << endl;
		cout << "File size: " << (double)fs::file_size(path) / (1024 * 1024) << " MB" << endl;
		cout << endl;

		// Compression //

		cout << "Compression in progress ... " << endl;

		if (Arithmetics == 1)
			Mesh->QuantizeCoordinates(MIO->GetQuantization());

		MIO->SetInput(Mesh);
		MIO->SetFileName(strout.c_str());	// output file name (instead of '.out')

		start = clock();
		MIO->Analyse();
		cout << "Analyse. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		start = clock();
		MIO->Synthetize();
		cout << "Synthetize. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		start = clock();
		MIO->Approximate();
		cout << "Approximate. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		start = clock();
		MIO->Write();
		cout << "Write. " << (double)(clock() - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

		Mesh->Delete();

		cout << "completed" << endl;
		///////////////// DO NOT NEED IF COMPRESSION END

		
		// Decompression //

		cout << "Decompression in progress ... " << endl;

		vtkMultiresolutionIO* MIO_d = vtkMultiresolutionIO::New();

		start = clock();
		MIO_d->SetInputFileName(basename.c_str()); // input file name setting
		MIO_d->SetFileName(strout.c_str());
		MIO_d->Read();
		end = clock();
		MIO_d->Delete();

		cout << "completed. " << (double)(end - start) / CLOCKS_PER_SEC << " sec. elapsed" << endl;

	}
}

