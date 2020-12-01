
#include <iostream>
#include <math.h>
#include <string>
//#include "add_func.cpp"


using namespace std;

//---------------------------------------------------------------------
//функции
void hello(int* ner_input_layer, int* n_hidd_layer, int* ner_hidd_layer, int* ner_out_layer,
			bool* a, double* b, double *lim_training,int *number_of_examples, bool* show_result_show,int* number_of_tests);
void wr_from_script(double* p, int n, string name_input, bool mod_files, int I_2);
void tr_fl(string n, bool* a);
void zero_1(double* p, int n);
void zero_2(double** p, int m, int n);
void random_1(double* p, int n);
void random_2(double** p, int m, int n);
void show_result(double* a,double *error, int n, int i1,bool ts, int* number_of_right);
void show(double* p, int n);
void wr_weights(double* p, int n, string name_1);
void wr_to_weights(double* a, string name_weights);

void one_side_way_ds(int* ner_hidd_layer, int* ner_input_layer, int* n_hidd_layer, int* ner_out_layer,
	double** ner_n_hidd_ds, double* ner_inp, double* ner_out, double* weight_input, double* weight_n_hidd,
	double* weight_output, long int* I1, long int* I2, long int* I3,
	int* length);

void back_error(int* ner_hidd_layer, int* ner_input_layer, int* n_hidd_layer, int* ner_out_layer,
	double** ner_n_hidd_ds, double* ner_inp, double* ner_out, double* weight_input, double* weight_n_hidd,
	double* weight_output, long int* I1, long int* I2, long int* I3,
	int* length, double* error_ner_out, double** error_ner_n_hidd_ds, double* result, long int* T1);

void change_weight(int* ner_input_layer, int* n_hidd_layer, int* ner_out_layer,
	double** ner_n_hidd_ds, double* ner_inp, double* weight_input, double* weight_n_hidd,
	double* weight_output, long int* I1, long int* I2, long int* I3,
	int* length, double* error_ner_out, double** error_ner_n_hidd_ds, double* normal_training);

string name_of_test_file(int test_number);
string name_of_training_file(int training_number);
void zero_ds(double** p, int elem, int layers);
int pr_mass(int* mass, int n);
double sigmoid(double x);
void show_result_ner(double* a, int n, int i1);

char latin(int num_l);
void db_show(double* a, int n, int h);
long int period(long  int a, int per);

void true_training_massive(bool* mass, int n);
//---------------------------------------------------------------------

//string name = "E:/files/c++/neiros/build-image-Desktop_Qt_5_15_1_MinGW_32_bit-Debug/photo.txt";
//string name_ad1 = "weights_1.txt";
//string name_ad2 = "weights_2.txt";


//нужна библиотека
//3 нейрона 3 скрытых слоя с 3 нейронами и двумя выходами



int main() {
	//доп переменные
	long int *I1 = new long int{}, *I2 = new long int{}, *I3 = new long int{};//для входа средней и выхода
	bool* training = new bool{};//тренировать
	int* training_number = new int{}; //training_number == number_of_examples?
	int* number_of_examples = new int{};
	int* number_of_tests = new int{};
	long int *T1 = new long int{};//для поседовательного считавания нужной части результатов
	bool ts{true};//делает так, чтобы в втором главном вложенном массиве четко T1 переходило до номера считвания следующих элементов(неактуально)
	bool* show_result_show = new bool{true};

	int* number_of_right = new int{};

	bool* training_final = new bool{ true };
	int* training_example = new int{};
	int* test_training_examples = new int{};



	double* lim_training = new double{};
	long int* element_of_iterat = new long int{};

	//количества элементов массивов
	int* ner_input_layer = new int{}, *n_hidd_layer = new int{}, *ner_hidd_layer = new int{}, *ner_out_layer = new int{};
	double* normal_training = new double{};



	
	hello(ner_input_layer, n_hidd_layer, ner_hidd_layer, ner_out_layer, training, normal_training, lim_training,number_of_examples, show_result_show, number_of_tests);


	//массив длин каждого скрытого слоя
	int* length = new int[*n_hidd_layer];
	for (int i{}; i < *n_hidd_layer; i++) {
		if (i == 0) { length[i] = *ner_hidd_layer; }
		else { length[i] = length[i - 1] / 2; }
	}
	//

	//доп массивы
	bool* training_massive = new bool[*number_of_examples];
	true_training_massive(training_massive, *number_of_examples);

	//----------------------------
	cout << "layers: ";
	for (int i{}; i < *n_hidd_layer; i++) {
		cout << length[i] << " ";
	}
	cout << '\n';
	//---------------------------

	//массивы нейронов
	double* ner_inp = new double[*ner_input_layer], * ner_out = new double[*ner_out_layer];
	
	//new
	int a{0};//нужный счетчик для создания массива
	double **ner_n_hidd_ds = new double* [*ner_hidd_layer];
	for (int i{}; i < *ner_hidd_layer; i++) {
		ner_n_hidd_ds[i] = new double[*n_hidd_layer - a];
		if (i == (length[ (*n_hidd_layer - 1 ) - a] - 1)) { a++; }
	}
	
	//

	double* result = new double[(*ner_out_layer)*(*ner_out_layer)];

	//--------------------------------------------------------
	//массивы весов
	double* weight_input = new double[(*ner_input_layer) * length[0]], * weight_output = new double[length[*n_hidd_layer-1]*(*ner_out_layer)];
	double* weight_n_hidd = new double [pr_mass(length,*n_hidd_layer-1)];
	//ошибки каждого нейрона
	double* error_ner_inp = new double[*ner_input_layer], * error_ner_out = new double[*ner_out_layer];
	a = 0;// layers = *n_hidd_layer;
	double** error_ner_n_hidd_ds = new double* [*ner_hidd_layer];
	for (int i{}; i < *ner_hidd_layer; i++) {
		error_ner_n_hidd_ds[i] = new double[*n_hidd_layer - a];
		if (i == (length[ (*n_hidd_layer - 1 ) - a] - 1)) { a++; }
	}
	//
	zero_1(ner_out, *ner_out_layer);
	zero_ds(ner_n_hidd_ds, *ner_hidd_layer, *n_hidd_layer);//обнуление выхода, входа, и средней части
	zero_1(error_ner_inp, *ner_input_layer);
	zero_1(error_ner_out, *ner_out_layer);
	zero_ds(error_ner_n_hidd_ds, *ner_hidd_layer, *n_hidd_layer);
		
	//--------------------------------------------------------
	if (*training) {

		random_1(weight_input, (*ner_input_layer) * length[0]);//первые значения весов
		random_1(weight_output, (*ner_out_layer) * length[*n_hidd_layer - 1]);
		random_1(weight_n_hidd, pr_mass(length, *n_hidd_layer - 1));

		

		//show(ner_inp, *ner_input_layer);


		wr_from_script(result, (*ner_out_layer) * (*ner_out_layer), "results.txt",false,0);
		//show(result, (*ner_out_layer) * (*number_of_examples)); 
		
		//db_show(result, (*ner_out_layer) * (*ner_out_layer), 26);
		while (*training_final){
				cout << "\nEpoch = element done ";
				

				

				for (int i_i{}; i_i < *number_of_examples; i_i++) {
					if (*element_of_iterat == *number_of_examples) { *element_of_iterat = 0; }

					if (*T1 == ((*ner_out_layer) * (*ner_out_layer))) { *T1 = 0; }

					wr_from_script(ner_inp, *ner_input_layer, "material.txt", true, (i_i * (*ner_input_layer)));
					//db_show(ner_inp, *ner_input_layer, 40);



					

					 
					//прямой проход
					one_side_way_ds(ner_hidd_layer, ner_input_layer, n_hidd_layer, ner_out_layer, ner_n_hidd_ds,
						ner_inp, ner_out, weight_input, weight_n_hidd, weight_output, I1, I2, I3, length);
					//----------------------------------------------------------------------------------------------------
					
					if (ner_out[period(*element_of_iterat,*ner_out_layer) ] > * lim_training) {
						*test_training_examples = *test_training_examples + 1;//проверка на узнавание всех букв за один цикл

						if (training_massive[*element_of_iterat]) {
							cout << latin(period(*element_of_iterat, *ner_out_layer) );
							cout << '\a';
							training_massive[*element_of_iterat] = false;
							
						}

						if ( *test_training_examples == *number_of_examples){
							*training_final = false;
							break;
						}

					}

					if (*show_result_show) {
						show_result_ner(ner_out, *ner_out_layer, *element_of_iterat);
					}

					//-----------------------------------------------

					//обратный проход( вычисление ошибки)
					back_error(ner_hidd_layer, ner_input_layer, n_hidd_layer, ner_out_layer, ner_n_hidd_ds,
						ner_inp, ner_out, weight_input, weight_n_hidd, weight_output, I1, I2, I3, length, error_ner_out, error_ner_n_hidd_ds, 
						result, T1);

					*I1 = 0; *I2 = 0; *I3 = 0;

					//-----------------------------------------------


					//изменение значений весов
					change_weight( ner_input_layer, n_hidd_layer, ner_out_layer, ner_n_hidd_ds,
						ner_inp, weight_input, weight_n_hidd, weight_output, I1, I2, I3, length, error_ner_out, error_ner_n_hidd_ds,
						 normal_training);

					
					//-----------------------------------------------
					
					zero_1(error_ner_out, *ner_out_layer);
					zero_ds(error_ner_n_hidd_ds, *ner_hidd_layer, *n_hidd_layer);

					zero_1(ner_out, *ner_out_layer);
					zero_ds(ner_n_hidd_ds, *ner_hidd_layer, *n_hidd_layer);

					*I1 = 0; *I2 = 0; *I3 = 0;


					*T1 = *T1 + (*ner_out_layer);

					*element_of_iterat = *element_of_iterat + 1;
				}

				cout << '\t' <<*test_training_examples ;
				*test_training_examples = 0;

				if (*show_result_show) {
					cout << "==================================================================================\n";
					cout << "==================================================================================\n";
					cout << "==================================================================================\n";
				}

				
			}

		wr_weights(weight_input, (*ner_input_layer)* length[0],"weights_1.txt");
		
		wr_weights(weight_n_hidd, pr_mass(length, *n_hidd_layer - 1), "weights_2.txt");
		
		wr_weights(weight_output, length[*n_hidd_layer-1]*(*ner_out_layer),"weights_3.txt");

		cout << "\n-----------------------------------------------------------------------\n";
		cout << "updated weights saved in txt files: weights_1.txt  weights_2.txt  weights_3.txt\n";
		cout << "-----------------------------------------------------------------------\n";

	}
	else {
		//готовая - проверка

		wr_to_weights(weight_input, "weights_1.txt");
		if (*n_hidd_layer != 1) {
			wr_to_weights(weight_n_hidd, "weights_2.txt");
		}
		wr_to_weights(weight_output, "weights_3.txt");//ошибка в данной строчке

		for (int i_p{}; i_p < *number_of_tests; i_p++) {

			wr_from_script(ner_inp, *ner_input_layer, name_of_test_file(i_p), false, 0);

			one_side_way_ds(ner_hidd_layer, ner_input_layer, n_hidd_layer, ner_out_layer, ner_n_hidd_ds, ner_inp, ner_out,
				weight_input, weight_n_hidd, weight_output, I1, I2, I3, length);

			show_result(ner_out, error_ner_out, *ner_out_layer, i_p, false, number_of_right);



			*I1 = 0; *I2 = 0; *I3 = 0;
			zero_1(ner_out, *ner_out_layer);
			zero_ds(ner_n_hidd_ds, *ner_hidd_layer, *n_hidd_layer);

		}

		cout << "\nNUMBER OF RIGHT ANSWERS = " << *number_of_right << '\n';
	}

	
	system("pause");
	return 0;
}

