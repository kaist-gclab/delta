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

#include "vtkMultiresolutionIO.h"
#include <time.h>
#include <string>
#include <iostream>
#include <filesystem>

int main( int argc, char *argv[] )
{
	cout << "Wavemesh, a wavelet based 3D mesh compression"<<endl;
	cout << "Copyright (C) 2008 S. Valette @ CREATIS, France "<<endl;
    cout << "(http://www-creatis.insa-lyon.fr/~valette)"<<endl;
	cout << "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY."<<endl;
	cout << "This program uses the VTK Library (www.kitware.com) for visualization"<<endl;
	cout << "This program uses the Range encoder (www.compressconsult.com) for coding"<<endl<<endl;

	if(argc<2)
	{
		cout<<"Usage :"<<endl;
		cout<<"Compression : wavemesh c file [options]"<<endl;
		cout<<"Decompression : wavemesh d file [options]"<<endl;
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
		cout<<"For mor help on wavemesh, execute the program without any argument"<<endl;
	}

	namespace fs = std::filesystem;

	std::string folder = argv[2];
	for (const auto & entry : fs::directory_iterator(folder))
	{
		// for debug (seperate modified model)
		std::cout << entry.path() << std::endl;
		std::string strtemp = entry.path().string();
		strtemp = strtemp.substr(strtemp.size() - 5, strtemp.size());
		if (strtemp[0] == 'm')
			continue;

		vtkMultiresolutionIO *MIO = vtkMultiresolutionIO::New();

		// input file name setting
		std::string inputfile = entry.path().string();
		inputfile = inputfile.substr(0, inputfile.size() - 4);
		MIO->SetInputFileName(inputfile.c_str());

		int Arithmetics = 1;

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
					cout << "Lifting radius : " << Lifting << endl;
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
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
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


		if (strcmp(argv[1], "cmp") == 0)
		{
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outfile = MIO->GetInputFileName();
			outfile += ".out";
			MIO->SetFileName(outfile.c_str());

			clock_t start, end;

			//start = clock();
			MIO->Analyse();
			MIO->Synthetize();
			//end = clock();
			//printf("Analyse, Synthetize elapsed time: %lf\n", (double)(end - start));

			if (MIO->GetLifting())
			{
				// re-compression
				while (1)
				{
					MIO->RoiSetting();

					//start = clock();
					MIO->Approximate();
					//end = clock();
					//printf("Approximate elapsed time: %lf\n", (double)(end - start));

					// timing experiment
					//start = time(NULL);
					MIO->Write();
					//end = time(NULL);
					//printf("Write elapsed time: %lf\n", end - start);
					break;
				}
			}
			else
			{
				// re-compression
				while (1)
				{
					//start = clock();
					MIO->Approximate();
					//end = clock();
					//printf("Approximate elapsed time: %lf\n", (double)(end - start));

					MIO->RoiSetting();

					// timing experiment
					//start = clock();
					MIO->Write();
					//end = clock();
					//printf("Write elapsed time: %lf\n", (double)(end - start));
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
						cout << "Lifting radius : " << Lifting << endl;
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
		// Compress -> Modify it (from _m.obj file) -> Recompress -> Decompress
		if (strcmp(argv[1], "r") == 0)
		{
			// merge first (since it has no parent-child relations)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
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
					printf("Approximate elapsed time: %lf\n", (double)(end - start));

					// timing experiment
					start = time(NULL);
					MIO->Write();
					end = time(NULL);
					printf("Write elapsed time: %lf\n", end - start);
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
					printf("Write elapsed time: %lf\n", (double)(end - start));
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
						cout << "Lifting radius : " << Lifting << endl;
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

		// d : decompression -> mesh.ply ( WriteOutput = 1(all mesh) or 2(10 finest mesh) )
		if (strcmp(argv[1], "d") == 0)
		{
			std::string strout = entry.path().string();
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

		// test1 : vertex reordering -> obj file
		if (strcmp(argv[1], "test1") == 0)
		{
			// 1. load each file
			// 2. make mesh file (after rearranging vertex ids)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			/*std::string outfile = entry.path().string();
			outfile = outfile.substr(0, outfile.size() - 4);
			outfile += ".out";
			MIO->SetFileName(outfile.c_str());*/

			clock_t start, end;

			//start = clock();
			MIO->Analyse();
			MIO->Synthetize();
			MIO->ExportOBJ(MIO->GetFilter(0)->GetOutput());
			//end = clock();
			//printf("Analyse, Synthetize elapsed time: %lf\n", (double)(end - start));

			Mesh->Delete();
			MIO->Delete();
		}

		// test2 : (no lifting) reordered mesh + modified reordered mesh -> .out + info.txt
		if (strcmp(argv[1], "test2") == 0)
		{
			std::string strout = MIO->GetInputFileName();
			strout += ".out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			/*std::string outfile = entry.path().string();
			outfile = outfile.substr(0, outfile.size() - 4);
			outfile += ".out";
			MIO->SetFileName(outfile.c_str());*/

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			MIO->RoiSetting();

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test3 : (lifting) reordered mesh + modified reordered mesh -> .out + info.txt
		if (strcmp(argv[1], "test3") == 0)
		{
			// output file name (instead of '.out')
			std::string strout = MIO->GetInputFileName();
			strout += "_lifting.out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			MIO->RoiSetting();

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test4 : (lifting) reordered mesh + randomly modify(different size) -> .out + info.txt
		if (strcmp(argv[1], "test4") == 0)
		{
			// output file name (instead of '.out')
			std::string strout = MIO->GetInputFileName();
			strout += "_lifting_random.out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			for (int k=1; k<=1; k*=2)
				MIO->RoiRandom(k);

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test34 : (lifting) reordered mesh + randomly modify(different size) + other modified mesh -> .out + info.txt
		if (strcmp(argv[1], "test34") == 0)
		{
			// output file name (instead of '.out')
			std::string strout = MIO->GetInputFileName();
			strout += "_lifting_random.out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			for (int k = 1; k <= 1024; k *= 2)
				MIO->RoiRandom(k);

			MIO->RoiSetting();

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test5 : (lifting) reordered mesh + incresingly modify(different size) -> .out + info.txt
		if (strcmp(argv[1], "test5") == 0)
		{
			// output file name (instead of '.out')
			std::string strout = MIO->GetInputFileName();
			strout += "_lifting_increase.out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			MIO->RoiIncrease(100);

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test45 : (lifting) reordered mesh + randomly modify(different size) + incresingly modify(different size) -> .out + info.txt
		if (strcmp(argv[1], "test45") == 0)
		{
			// output file name (instead of '.out')
			std::string strout = MIO->GetInputFileName();
			strout += "_lifting_increase.out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			for (int k = 1; k <= 1024; k *= 2)
				MIO->RoiRandom(k);

			out.open(outinfo, ios::app);
			out << endl;
			out.close();

			MIO->RoiIncrease(50);

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test6 : files -> one info.txt (include filename, #vertices, #faces)
		if (strcmp(argv[1], "test6") == 0)
		{
			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());

			// input file name
			std::string inputfile = MIO->GetInputFileName();

			std::string outinfo = entry.path().parent_path().string() + "\\all_infos.txt";
			std::ofstream out;

			out.open(outinfo, ios::app);
			out << inputfile << "  " << Mesh->GetNumberOfPoints() << "  " << Mesh->GetNumberOfCells() << endl;
			out.close();

			Mesh->Delete();
		}

		// test7 : (lifting) reordered mesh + k% local modification (start at random initial point, within [k,k+0.5]% ROI -> .out + info.txt
		if (strcmp(argv[1], "test7") == 0)
		{
			// output file name (instead of '.out')
			std::string strout = MIO->GetInputFileName();
			strout += "_lifting_local.out";
			MIO->SetFileName(strout.c_str());

			// 1. load each file
			// 2. compression
			// 3. make output
			// 4. check times (Approximate time, Write time, number of modified vertices)
			vtkSurface *Mesh = vtkSurface::New();
			cout << "Load : " << entry.path() << endl;
			Mesh->CreateFromFile(entry.path().string().c_str());
			Mesh->DisplayMeshProperties();

			if (Arithmetics == 1)
				Mesh->QuantizeCoordinates(MIO->GetQuantization());

			MIO->SetInput(Mesh);

			std::string outinfo = MIO->GetInputFileName();
			outinfo += "_info.txt";
			std::ofstream out;

			out.open(outinfo, ios::trunc);
			out.close();

			clock_t start, end;

			start = clock();
			MIO->Analyse();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Synthetize();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			start = clock();
			MIO->Approximate();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			MIO->RoiLocal(5);

			start = clock();
			MIO->Write();
			end = clock();

			out.open(outinfo, ios::app);
			out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
			out.close();

			Mesh->Delete();
			MIO->Delete();
		}

		// test8 : (various lifting radius k) reordered mesh + full modify -> .out + info.txt
		if (strcmp(argv[1], "test8") == 0)
		{
			int disks[] = { 0,1,2,3 };
			for (int r = 0; r < sizeof(disks)/sizeof(int); r++)
			{
				vtkMultiresolutionIO* MIO = vtkMultiresolutionIO::New();

				// input file name setting
				std::string inputfile = entry.path().string();
				inputfile = inputfile.substr(0, inputfile.size() - 4);
				MIO->SetInputFileName(inputfile.c_str());

				int Arithmetics = 1;

				///// Set Lifting Radius /////
				MIO->SetLifting(1);
				MIO->SetLiftingRadius(disks[r]);
				char radius[2];
				sprintf(radius, "%d", disks[r]);

				// output file name (instead of '.out')
				std::string strout = MIO->GetInputFileName();
				strout += "_lifting_random_";
				strout += radius;	// file name + radius
				strout += ".out";
				MIO->SetFileName(strout.c_str());

				// 1. load each file
				// 2. compression
				// 3. make output
				// 4. check times (Approximate time, Write time, number of modified vertices)
				vtkSurface* Mesh = vtkSurface::New();
				cout << "Load : " << entry.path() << endl;
				Mesh->CreateFromFile(entry.path().string().c_str());
				Mesh->DisplayMeshProperties();

				if (Arithmetics == 1)
					Mesh->QuantizeCoordinates(MIO->GetQuantization());

				MIO->SetInput(Mesh);

				std::string outinfo = MIO->GetInputFileName();
				outinfo += "_info_";
				outinfo += radius;	// file name + radius
				outinfo += ".txt";
				std::ofstream out;

				out.open(outinfo, ios::trunc);
				out.close();

				clock_t start, end;

				start = clock();
				MIO->Analyse();
				end = clock();

				out.open(outinfo, ios::app);
				out << "Analyse elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
				out.close();

				start = clock();
				MIO->Synthetize();
				end = clock();

				out.open(outinfo, ios::app);
				out << "Synthetize elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
				out.close();

				start = clock();
				MIO->Approximate();
				end = clock();

				out.open(outinfo, ios::app);
				out << "Approximate elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
				out.close();

				MIO->RoiRandom(1);

				start = clock();
				MIO->Write();
				end = clock();

				out.open(outinfo, ios::app);
				out << "Write elapsed time: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
				out.close();

				Mesh->Delete();
				MIO->Delete();
			}
			
		}
	}
}

