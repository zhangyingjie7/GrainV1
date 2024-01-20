#include "GrainV1_Gamma.h"

int GrainV1_New()
{
	ofstream ofile("GrainV1_New_Gamma.txt");
	double H_corr_value[2][16] = { {0, 0, 0, 0, 0, -8, 0, 8 ,0, 8, 0, -8, -8, 8, -8, 8},\
	{0, -8, 0, 8, -8, -8, -8, -8, 0, 0, 0, 0, 0, -8, 0, 8} };
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			H_corr_value[i][j] = H_corr_value[i][j] / 32;
		}
	}
	/*0,9,21,37,52,60,62,80*/
	int SpanW[32][8] = { {0, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 0, 0, 0, 0, 1, 2} ,{0, 0, 0, 2, 0, 0, 0, 8} ,{0, 0, 0, 2, 0, 0, 1, 10},\
	{0, 0, 1, 0, 0, 4, 0, 0} ,{0, 0, 1, 0, 0, 4, 1, 2} ,{0, 0, 1, 2, 0, 4, 0, 8} ,{0, 0, 1, 2, 0, 4, 1, 10} ,\
	{0, 2, 0, 0, 8, 0, 0, 0} ,{0, 2, 0, 0, 8, 0, 1, 2} ,{0, 2, 0, 2, 8, 0, 0, 8} ,{0, 2, 0, 2, 8, 0, 1, 10} ,\
	{0, 2, 1, 0, 8, 4, 0, 0} ,{0, 2, 1, 0, 8, 4, 1, 2} ,{0, 2, 1, 2, 8, 4, 0, 8} ,{0, 2, 1, 2, 8, 4, 1, 10} ,\
	{2, 0, 4, 0, 0, 0, 0, 0} ,{2, 0, 4, 0, 0, 0, 1, 2} ,{2, 0, 4, 2, 0, 0, 0, 8} ,{2, 0, 4, 2, 0, 0, 1, 10} ,\
	{2, 0, 5, 0, 0, 4, 0, 0} ,{2, 0, 5, 0, 0, 4, 1, 2} ,{2, 0, 5, 2, 0, 4, 0, 8} ,{2, 0, 5, 2, 0, 4, 1, 10} ,\
	{2, 2, 4, 0, 8, 0, 0, 0} ,{2, 2, 4, 0, 8, 0, 1, 2} ,{2, 2, 4, 2, 8, 0, 0, 8} ,{2, 2, 4, 2, 8, 0, 1, 10} ,\
	{2, 2, 5, 0, 8, 4, 0, 0} ,{2, 2, 5, 0, 8, 4, 1, 2} ,{2, 2, 5, 2, 8, 4, 0, 8} ,{2, 2, 5, 2, 8, 4, 1, 10} };

	int G_corr_list[16][12] = { {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0} ,\
	{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0} ,{0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0} ,\
	{0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0} ,{0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0} ,\
	{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0} ,{0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0} ,\
	{0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0} ,{0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0} ,\
	{0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0} ,{0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0} };

	double G_corr_Value[16] = { 5.23146559316956e-06, 2.0613540527847363e-06, 4.993274842490791e-06, 1.971279061763198e-06,\
		5.161289209354436e-07, -5.894669357076054e-07, 4.918340437143343e-07, -5.628255621559219e-07,\
		- 8.160088782460662e-07, -5.892622993997065e-07, -7.677097073610639e-07, -5.630301984638209e-07, \
		- 5.159242846275447e-07, 5.892622993997065e-07, -4.920386800222332e-07, 5.630301984638209e-07 };

	//double threshold = pow(2, -36);
	double threshold = pow(2, -35.8448);
	int count = 0;
	double corr;
	double corr_temp;
	double minimum = 1;
	double maximum = 0;
	double got;

	//DWORD start_time = GetTickCount();
	int i14_array[3] = { 9,7,5 };
	int i28_array[4] = { 9,7,5,1 };
	int i45_array[3] = { 9,7,5 };


	for (int i0 = 0; i0 < 8; i0++)
	{
		cout << "i0=" << i0 << endl;
		int tmp_i0 = (i0 & 1) | ((i0 >> 1 & 1) << 2) | ((i0 >> 2 & 1) << 3);
		for (int i9 = 0; i9 < 8; i9++)
		{
			int tmp_i9 = (i9 & 1) | ((i9 >> 1 & 1) << 2) | ((i9 >> 2 & 1) << 3);
			for (int &i14 :i14_array)
			{
				for (int i21 = 0; i21 < 4; i21++)
				{
					int tmp_i21 = ((i21 & 1) << 1) | ((i21 >> 1 & 1) << 3);
					for (int &i28:i28_array)
					{
						for (int i33 = 15; i33 < 16; i33++)
						{
							for (int i37 = 0; i37 < 8; i37++)
							{
								int tmp_i37 = (i37 & 1) | ((i37 >> 1 & 1) << 2) | ((i37 >> 2 & 1) << 3);
								for (int &i45:i45_array)
								{
									for (int i55 = 0; i55 < 8; i55++)
									{
										int tmp_i55 = i55;
										for (int i60 = 0; i60 < 8; i60++)
										{
											int tmp_i60 = (i60 & 1) | ((i60 >> 1 & 1) << 1) | ((i60 >> 2 & 1) << 3);
											for (int i62 = 0; i62 < 8; i62++)
											{
												int tmp_i62 = ((i62 & 1) << 1) | ((i62 >> 1 & 1) << 2) | ((i62 >> 2 & 1) << 3);
												for (int i80 = 0; i80 < 4; i80++)
												{
													int tmp_i80 = (i80 & 1) | ((i80 >> 1 & 1) << 2);
													//cout << i14<<" "<<i28<<" "<<i45 << endl;
													int linear_mask[12] = { tmp_i0, tmp_i9,i14,tmp_i21,i28 ,i33,tmp_i37 ,i45 ,tmp_i55,tmp_i60, tmp_i62 ,tmp_i80 };
													corr = 0;
													for (int(*gamma_tz)[12] = G_corr_list; gamma_tz != G_corr_list + 16; gamma_tz++)//v
													{
														for (int(*d_w)[8] = SpanW; d_w != SpanW + 32; d_w++) // w
														{
															int *p2 = *d_w;
															corr_temp = G_corr_Value[gamma_tz - G_corr_list];

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][0]][linear_mask[0] ^ (*p2)]; /*i0*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][1]][linear_mask[1] ^ (*p2 + 1)]; /*i9*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][2]][linear_mask[2]]; /*i14*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][3]][linear_mask[3] ^ (*p2 + 2)]; /*i21*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][4]][linear_mask[4]]; /*i28*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][5]][linear_mask[5]]; /*i33*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][6]][linear_mask[6] ^ (*p2 + 3)]; /*i37*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][7]][linear_mask[7]]; /*i45*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][8]][linear_mask[8] ^ (*p2 + 4)]; /*i52*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][9]][linear_mask[9] ^ (*p2 + 5)]; /*i60*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][10]][linear_mask[10] ^ (*p2 + 6)]; /*i62*/
															corr_temp = corr_temp * got;

															got = H_corr_value[G_corr_list[gamma_tz - G_corr_list][11]][linear_mask[11] ^ (*p2 + 7)]; /*i80*/
															corr_temp = corr_temp * got;

															corr = corr + corr_temp;
														}
													}

													//if (corr != 0)
													//{
													//	if (abs(corr) < abs(minimum))
													//	{
													//		minimum = corr;
													//		cout << "MIN " << minimum / abs(minimum) * 2 << "^" << log2(abs(minimum)) << endl;
													//	}
													//	//cout << corr / abs(corr) * 2 << "^" << log2(abs(corr)) << endl;

													//	if (abs(corr) > abs(maximum))
													//	{
													//		maximum = corr;
													//		cout << "MAX " << maximum / abs(maximum) * 2 << "^" << log2(abs(maximum)) <<" (" <<i14 << " " << i28 << " " << i45 <<")"<< endl;
													//	}
													//}

													if (abs(corr) >= threshold)
													{
														cout<<" (" << i14 << " " << i28 << " " << i45 << ")" << endl;
														ofile << "(";
														for (int cf = 0; cf < 11; cf++)
														{
															ofile << linear_mask[cf] << ",";
														}
														ofile << linear_mask[11] << "),  " << corr / abs(corr) * 2 << "^" << log2(abs(corr)) << "  (" << i14 << " " << i28 << " " << i45 << ")" << endl;
														//ofile << linear_mask[11] << "),  " << corr << endl;
														count = count + 1;
													}
												}
											}
										}
									}
								}

							}
						}
					}
				}
			}
		}
	}
	ofile.close();
	//DWORD end_time = GetTickCount();
	//cout << "The run time is:" << (end_time - start_time) << "ms!" << endl;
	//cout << "Min correlation is: " << minimum / abs(minimum) * 2 << "^" << log2(abs(minimum)) << endl;
	//cout << "Max correlation is: " << maximum / abs(maximum) * 2 << "^" << log2(abs(maximum)) << endl;
	cout << "There are " << count << " correlations with abs >= " << 2 << "^" << log2(threshold) << endl;
	//cout << "There are " << count << " correlations with abs > " << threshold << endl;
	cout << "count = 2^" << (log2(abs(count))) << endl;
	//system("pause");
	return 0;
}
