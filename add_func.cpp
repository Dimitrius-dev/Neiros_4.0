#include <string>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

char latin(int num_l) {
	return num_l + 65;
}


void show_result_ner(double* a, int n, int i1) {
	cout << "\n--------------------------- => ";
	cout << "out number = " << i1 + 1 << '\n';
	for (int i{}; i < n; i++) {
		cout << "\nresult output neuron nam." << latin(i) << " = " << a[i];
	}
	cout << "\n";
}



void show_result(double* a, double* error, int n, int i1, bool ts,int *number_of_right) {//i1 счётчик
	double mx{}; int* i_mx = new int{};
	cout << "---------------------------------------------\n";
	cout << "out number = " << i1 + 1 << '\n';
		if (ts) {
			for (int i{}; i < n; i++) {
				cout << "\tresult output neuron nam." << i + 1 << " = " << a[i] << "        error = " << error[i] << '\n';
			}
		}
		else {
			for (int i{}; i < n; i++) {
				cout << "\tresult output neuron nam." << latin(i) << " number" << i << " = " << a[i] << '\n';//i + 1
			}

			for (int i{}; i < n; i++) {
				if (a[i] > mx) { mx = a[i]; *i_mx = i; }
			}

			cout << "\nmachine think that this number is =  " << latin (*i_mx) << '\n';
			if (latin(*i_mx) == latin(i1)) {
				cout << "\n=========== RIGHT ===========\n";
				*number_of_right = *number_of_right + 1;
			}
			else {
				cout << "\n=========== FALSE ===========\n";
			}
		}

	//нахождение максимального


}



void zero_1(double* p, int n) {
	for (int i{}; i < n; i++) { p[i] = 0; }
}

void zero_2(double** p, int m, int n) {
	for (int i{}; i < m; i++) {
		for (int j{}; j < n; j++) {
			p[i][j] = 0.;
		}
	}
}

void zero_ds(double** p, int elem,int layers) {
	int elem_1 = elem;

	for (int i{}; i < layers; i++) {
		for (int j{}; j < elem_1; j++) {
			p[j][i] = 0.;
		}
		elem_1 /= 2;
	}
}

void random_1(double* p, int n) {
	for (int i{}; i < n; i++) {
		p[i] = ((rand() % 10) / 10.) - 0.5;
	}
}


//old
void random_2(double** p, int m, int n) {
	for (int i{}; i < m; i++) {
		for (int j{}; j < n; j++) {
			p[i][j] = ((rand() % 10) / 10.) - 0.5;
		}
	}
}
//
//new
void random_ds(double** p, int m, int n) {
	int elem_1 = m;//m - down on massive n - rigth on massive
	for (int i{}; i < n; i++) {
		for (int j{}; j < elem_1; j++) {
			p[j][i] = ((rand() % 10) / 10.) - 0.5;
		}
		elem_1 /= 2;
	}
}
//
void show(double* p, int n) {
	for (int i{}; i < n; i++) {
		cout << p[i] << '\n';
	}
}

void show_db_2(double** p, int m, int n) {
	for (int i{}; i < m; i++) {
		for (int j{}; j < n; j++) {
			cout << p[j][i] << '\n';
		}
	}
}

double sigmoid(double x) {
	return (1 / (1 + exp(-x)));
}

void wr_from_script(double* p, int n, string name_input, bool mod_files, int I_2) {
	ifstream output; long int I_1{ 0 }, I_ch{}; char buf{};
	output.open(name_input);
	if (output.is_open()) {

		if (!mod_files) {
			while (output.get(buf)) {
				if (I_1 == n) { break; }
				if (buf == '\n') { continue; }
				p[I_1] = buf - 48.; I_1++;
			}
		}
		else {
			while (output.get(buf)) {
				if (buf == '\n') { continue; }
				if (I_ch == I_2) {
					if (I_1 == n) { break; }
					p[I_1] = buf - 48.; I_1++;
				}
				else {
					I_ch++;
				}

			}
		}
	}
	else { cout << "\nERROR OPEN FILE_1\n"; }
	output.close();
}


void wr_to_weights(double* a, string name_weights) {
	ifstream in_wh; long int I{}; string buffer = "";
	in_wh.open(name_weights);
	if (in_wh.is_open()) {
		while (!in_wh.eof()) {
			in_wh >> buffer;
			a[I] = stod(buffer); I++;
		}
	}
	else { cout << "\nERROR OPEN FILE_2\n"; }
	in_wh.close();
}



void wr_weights(double* p, int n, string name_1) {
	ofstream file; char buf{}; int m{};
	file.open(name_1);
	if (file.is_open()) {
		for (int i{}; i < n; i++) {
			if (m == 50) {
				file << p[i] << '\n';
				m = 0;
			}
			else {
				file << p[i] << " ";
			}

		}
	}
	else { cout << "\nERROR OPEN FILE_3\n"; }
	file.close();
}

void tr_fl(string n, bool* a) {
	if ((n == "Yes") || (n == "yes") || (n == "YES") || (n == "y") || (n == "Y")) { *a = true; }
	if ((n == "No") || (n == "no") || (n == "NO") || (n == "n") || (n == "N")) { *a = false; }
}



int pr_mass(int *mass,int n) {
	int mx{};
	for (int i{}; i < n; i++) {
		mx += (mass[i] * mass[i + 1]);
	}
	return mx;
	
}

void one_side_way_ds(int* ner_hidd_layer, int* ner_input_layer, int* n_hidd_layer, int* ner_out_layer,
	double** ner_n_hidd_ds, double* ner_inp, double* ner_out, double* weight_input, double* weight_n_hidd,
	double* weight_output, long int* I1, long int* I2, long int* I3,
	int *length) {

	for (int i1{}; i1 < length[0]; i1++) {//вход
		for (int j1{}; j1 < *ner_input_layer; j1++) {
			ner_n_hidd_ds[i1][0] += (ner_inp[j1] * weight_input[*I1]);// cout << weight_input[I1]<<'\n';
			*I1 = *I1 + 1;
		}
		ner_n_hidd_ds[i1][0] = sigmoid(ner_n_hidd_ds[i1][0]);// cout << ner_n_hidd[i1][0] << '\n';
	}

	//cout << '\n';

	for (int i1{ 1 }; i1 < *n_hidd_layer; i1++) {//среднее //(int i1{ 1 }; i1 <= *n_hidd_layer; i1++)
		for (int j1{}; j1 < length[i1]; j1++) {
			for (int k1{}; k1 < length[i1-1]; k1++) {
				ner_n_hidd_ds[j1][i1] += (ner_n_hidd_ds[k1][i1 - 1] * weight_n_hidd[*I2]); //cout << weight_n_hidd[I2] << '\n';
				*I2 = *I2 + 1;
			}
			ner_n_hidd_ds[j1][i1] = sigmoid(ner_n_hidd_ds[j1][i1]); //cout << ner_n_hidd[j1][i1] << '\n';
		}
	}

	//cout << '\n';

	for (int i1{}; i1 < *ner_out_layer; i1++) {//выход
		for (int j1{}; j1 < length[*n_hidd_layer-1]; j1++) {
			ner_out[i1] += (ner_n_hidd_ds[j1][(*n_hidd_layer) - 1] * weight_output[*I3]); ///cout << weight_output[I3] << '\n'; 
			*I3 = *I3 + 1;
		}
		ner_out[i1] = sigmoid(ner_out[i1]);// cout << ner_out[i1] << '\n';
	}

}


void back_error(int* ner_hidd_layer, int* ner_input_layer, int* n_hidd_layer, int* ner_out_layer,
	double** ner_n_hidd_ds, double* ner_inp, double* ner_out, double* weight_input, double* weight_n_hidd,
	double* weight_output, long int* I1, long int* I2, long int* I3,
	int* length, double * error_ner_out,double ** error_ner_n_hidd_ds,double * result, long int *T1) {

	for (int i1{}; i1 < (*ner_out_layer); i1++) {
		error_ner_out[i1] = (result[i1 + *T1] - ner_out[i1]) * (ner_out[i1]) * (1 - ner_out[i1]);

	}


	//распространение ошибки
	for (int i1{ (*ner_out_layer) - 1 }; i1 >= 0; i1--) {//выход
		for (int j1{ (length[*n_hidd_layer - 1]) - 1 }; j1 >= 0; j1--) {
			*I3 = *I3 - 1;
			error_ner_n_hidd_ds[j1][*n_hidd_layer - 1] += error_ner_out[i1] * weight_output[*I3] * (ner_n_hidd_ds[j1][(*n_hidd_layer) - 1]) * (1 - (ner_n_hidd_ds[j1][(*n_hidd_layer) - 1]));
		}
	}

	for (int i1{ *n_hidd_layer - 1 }; i1 >= 1; i1--) {//среднее
		for (int j1{ length[i1] - 1 }; j1 >= 0; j1--) {
			for (int k1{ length[i1 - 1] - 1 }; k1 >= 0; k1--) {
				*I2 = *I2 - 1;
				error_ner_n_hidd_ds[k1][i1 - 1] += error_ner_n_hidd_ds[j1][i1] * weight_n_hidd[*I2] * (ner_n_hidd_ds[k1][i1 - 1]) * (1 - (ner_n_hidd_ds[k1][i1 - 1]));//ner_n_hidd[j1][i1] += (ner_n_hidd[k1][i1 - 1] * weight_n_hidd[I2]); I2++;
			}
		}
	}
}

void change_weight( int* ner_input_layer, int* n_hidd_layer, int* ner_out_layer,
	double** ner_n_hidd_ds, double* ner_inp, double* weight_input, double* weight_n_hidd,
	double* weight_output, long int* I1, long int* I2, long int* I3,
	int* length, double* error_ner_out, double** error_ner_n_hidd_ds,double* normal_training) {
	for (int i1{}; i1 < length[0]; i1++) {//вход
		for (int j1{}; j1 < *ner_input_layer; j1++) {
			weight_input[*I1] += ner_inp[j1] * error_ner_n_hidd_ds[i1][0] * (*normal_training); *I1 = *I1 + 1; //ner_n_hidd[i1][0] += (ner_inp[j1] * weight_input[I1]); I1++;
		}

	}

	for (int i1{ 1 }; i1 < *n_hidd_layer; i1++) {//среднее //(int i1{ 1 }; i1 <= *n_hidd_layer; i1++)
		for (int j1{}; j1 < length[i1]; j1++) {
			for (int k1{}; k1 < length[i1 - 1]; k1++) {
				weight_n_hidd[*I2] += ner_n_hidd_ds[k1][i1 - 1] * error_ner_n_hidd_ds[j1][i1] * (*normal_training); *I2 = *I2 + 1;//ner_n_hidd[j1][i1] += (ner_n_hidd[k1][i1 - 1] * weight_n_hidd[I2]); I2++;
			}

		}
	}

	//cout << '\n' << length[*n_hidd_layer - 1] << '\n';

	for (int i1{}; i1 < *ner_out_layer; i1++) {//выход
		for (int j1{}; j1 < length[*n_hidd_layer - 1]; j1++) {
			weight_output[*I3] += ner_n_hidd_ds[j1][(*n_hidd_layer) - 1] * error_ner_out[i1] * (*normal_training); *I3 = *I3 + 1;
		}

	}
}


void hello(int* ner_input_layer, int* n_hidd_layer, int* ner_hidd_layer, int* ner_out_layer,
	bool* a, double* b,double  *lim_training, int* number_of_examples, bool* show_result_show, int* number_of_tests) {
	string tr = "";
	cout << "--------------------------------------------------------------------\n";
	cout << "-----------------------------H-E-L-L-O------------------------------\n";
	cout << "--------------------------------------------------------------------\n";
	cout << "---------------------This-is-neuronet-by-Dmitrii-B.-----------------\n";
	cout << "--------------------------------------------------------------------\n";
	cout << "-----------------------------version-3.1----------------------------\n";
	cout << "--------------------------------------------------------------------\n";
	cout << "--------------------------------------------------------------------\n";

	cout << "Enter n  neurons of input_layer = "; cin >> *ner_input_layer;
	cout << "Enter n for hidden_layers = "; cin >> *n_hidd_layer;
	cout << "Enter n neurons of hidden_layer = "; cin >> *ner_hidd_layer;
	cout << "Enter n neurons of output_layer = "; cin >> *ner_out_layer;
	cout << "Enter number of tests = "; cin >> *number_of_tests;

	cout << "Enter training mode(yes/no) = "; cin >> tr; tr_fl(tr, a);
	if (*a == true) {
		cout << "==============================================================\n";
		cout << "   Enter normal training(0.0 - 1.0) = "; cin >> *b;
		cout << "   Enter lim training(7.0 - 0.99 examle) = "; cin >> *lim_training;
		cout << "   Enter set of number of examples = "; cin >> *number_of_examples; *number_of_examples = (*number_of_examples)*(*ner_out_layer);
		cout << "   Show process information - slow working(yes/no) = "; cin >> tr; tr_fl(tr, show_result_show);
	}
}

string name_of_training_file(int training_number) {
	return ( to_string(training_number + 1) + ".txt" );
}

string name_of_test_file(int test_number) {
	return ("test_" + to_string(test_number + 1) + ".txt");
}

	
void db_show(double* a,int n,int h) {
	long int I{0};
	for (int i{}; i < n; i++) {
		if (I == h) { cout << '\n'; I = 0; }
		cout << a[i];
		I++;
		
	}
	cout << '\n';
	cout << '\n';
}


/*
	if(*show_result_show) {
		show_result(ner_out, error_ner_out, *ner_out_layer, i, true, number_of_right);
	}
*/

void true_training_massive(bool *mass,  int n) {
	for (int i{}; i < n; i++) { 
		mass[i] = true;
	}
}

long int period(long  int a, int per) {
	if (a >= per) {
		return a % per;
	}
	else {
		return a;
	}
}