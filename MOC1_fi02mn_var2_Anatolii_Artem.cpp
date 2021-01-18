#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <numeric>


void collect_the_csv_values(std::vector<float>& pM, std::vector<float>& pK, std::vector<float>& pC, std::vector<std::vector<uint32_t>>& table)
{
	std::fstream csv_read_probs("/home/ftl/moc1/vars/prob_06.csv", std::ios::in);
	std::fstream csv_read_table("/home/ftl/moc1/vars/table_06.csv", std::ios::in);

	std::string temp;
////////////////////////////reading prob_XX.csv ////////////////////////
	getline(csv_read_probs, temp);
	std::cout << temp << std::endl;
	for(size_t i = 0; i < 20; ++i)
	{
		pM.push_back(std::stof(temp));
		if(i != 19)
		{
			temp.erase(temp.begin(), temp.begin()+5);
		}
	}


	getline(csv_read_probs, temp);
	std::cout << temp << std::endl;
	for(size_t i = 0; i < 20; ++i)
	{
		pK.push_back(std::stof(temp));
		if(i != 19)
		{
			temp.erase(temp.begin(), temp.begin()+5);
		}
	}


//////////////////////////////////reading table_XX ////////////////////////////////
	for(size_t i = 0; i < 20; ++i)
	{
		getline(csv_read_table, temp);
		for(size_t j = 0; j < 20; ++j)
		{
			table[i][j] = (std::stoul(temp));
			if(j != 19)
			{
				temp.erase(0 , temp.find_first_of(',')+1);
			}

		}
	}

	for(auto it = table.begin(); it != table.end(); ++it)
	{
		for(auto it1 = (*it).begin(); it1 != (*it).end(); ++it1)
		{
			std::cout << *it1 << ", ";
		}
		std::cout << std::endl;
	}

///////////////////////////////////calculation of std::vector<float> pC /////////////////////////////

	std::vector<std::vector<float>> pMandK(20, std::vector<float>(20,0));
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			pMandK[i][j] = (pM[i] * pK[j]);
		}
	}
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			pC[table[i][j]] += pMandK[j][i]; 
		}
	}

	


	csv_read_probs.close();
	csv_read_table.close();


}

void create_pMifC_csv(std::vector<float> pM, std::vector<float> pK, std::vector<float> pC, std::vector<std::vector<uint32_t>> table, std::vector<std::vector<float>>& pMandC, std::vector<std::vector<float>>& pMifC)
{
	std::fstream csv_write_pMifC("/home/ftl/moc1/pMifC.csv", std::fstream::out);

//////////////////////////////calculation of pMandC /////////////////////////////

	std::vector<std::vector<float>> pMandK(20, std::vector<float>(20,0));
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			pMandK[i][j] = (pM[i] * pK[j]);
		}
	}
	
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			pMandC[i][table[j][i]] += pMandK[i][j]; 
		}
	}
	
//////////////////////////////calculation of pMifC////////////////////////////

	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			pMifC[i][j] = (pMandC[i][j] / pC[j]);
		}
	}

	auto copypMifC = pMifC;
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			copypMifC[j][i] = pMifC[i][j];
			std::cout << pMifC[i][j] << ", ";
		}
		std::cout << std::endl;
	}	

	
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			csv_write_pMifC << copypMifC[i][j] << ',';
		}
		csv_write_pMifC << '\n';
	}

	csv_write_pMifC.close();
}

std::vector<std::pair<uint32_t,uint32_t>> Deterministic(std::vector<std::vector<uint32_t>> table, std::vector<std::vector<float>> pMandC, std::vector<std::vector<float>> pMifC)
{
	std::vector<std::pair<uint32_t,uint32_t>> CT_and_PT(20, std::pair<uint32_t, uint32_t>(0,0));
	float vitrati = 0.0;
	uint32_t answ = 0;
	/* //DEBUG
	
	std::cout << std::endl << "------------------------------------------------" << std::endl;
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			std::cout << pMandC[i][j] << ",  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << "------------------------------------------------" << std::endl;
	std::cout << std::endl << "------------------------------------------------" << std::endl;
	for(size_t i = 0; i < 20; ++i)
	{
		for(size_t j = 0; j < 20; ++j)
		{
			std::cout << pMifC[i][j] << ",  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << "------------------------------------------------" << std::endl;
	*/

	for(size_t col = 0; col < table.size(); ++col)
	{
		answ = 0;
		for(size_t row = 0; row < table.size(); ++row)
		{
			if(pMifC[answ][col] < pMifC[row][col])
			{
				answ = row;//std::max_element if transponate the matrix
			}

		}
		
		if(answ != table[col][answ])
		{
			vitrati = vitrati + pMandC[answ][col];
		}
		
		CT_and_PT[col] = {col, answ};
	}
	std::cout << "VIT " <<  vitrati << std::endl; 
	std::cout << "1.0 - VIT " << 1.0 - vitrati << std::endl; 
	return CT_and_PT;
}

std::vector<std::pair<uint32_t, uint32_t>> Stohastic(std::vector<std::vector<uint32_t>> table, std::vector<std::vector<float>> pMandC, std::vector<std::vector<float>> pMifC)
{
	std::vector<std::pair<uint32_t, uint32_t>> ret(20, std::pair<uint32_t, uint32_t>(0,0));
	float vitrati = 0.0;
	uint32_t answ = 0;
	float L = 0.0;
		float cos = 0.0;
		std::vector<std::vector<float>> L_matrix(20,std::vector<float>(20,1.0));
		std::vector<std::vector<uint32_t>> Answ_coord (20, std::vector<uint32_t>(0,0));
	for(size_t col = 0; col < table.size(); ++col)
	{
		ret[col].first = col;
		answ = 0;
	//	vitrati = 0.0;
		uint32_t amount = 0;
		auto biggest = pMifC[0][col];
		uint32_t it = 0;
		for(size_t row = 0; row < table.size(); ++row)
		{
			if(pMifC[row][col] > pMifC[answ][col])
			{
				answ = row;
				amount = 1;
				biggest = pMifC[row][col]; 
				it = row;
			}
				if(pMifC[row][col] == pMifC[answ][col])
				{
					amount += 1;
				}
		}

		std::cout << "AMOUNT " <<amount <<  std::endl;
		//std::cout << "BIGGEST " << biggest <<  std::endl;

		for(size_t row = 0; row < table.size(); ++row)
		{

			if(pMifC[row][col] == pMifC[it][col])
			{
				Answ_coord[col].push_back(row); 
			}
			
		}

	
/*

		for(auto it = Answ_coord.begin(); it != Answ_coord.end(); ++it)
		{
			for(auto iter = (*it).begin(); iter != (*it).end(); ++iter)
			{
				std::cout << *iter << ", ";
			}
			std::cout << std::endl;
		}

*/

		L = 0.0;
		
	//	L = (20-amount) *  (((float)(1.0) / (float)(amount)));
	// L matrix 20x20
	// L x pMandC ; sum (LxP; L)
		vitrati +=   ((L) * ( pMandC[it][col]));

		ret[col].second = it;



		
//	std::cout << "prob for " << col << " is " << ((float)(1.0) / ((float)amount)) << std::endl;
//	std::cout << "vitrati for " << col << " is " << vitrati << std::endl;



	}	
		for(size_t i = 0; i < table.size(); ++i)
		{
			for(size_t j = 0; j < Answ_coord[i].size(); ++j)
			{
				L_matrix[i][Answ_coord[i][j]] = 1.0 - ( 1.0 / (float)(Answ_coord[i].size() ) ); 
				//std::cout << Answ_coord[i][j] << ", ";

			}
		}
	
		for(size_t i = 0; i < table.size(); ++i)
		{
			for(size_t j = 0; j < Answ_coord[i].size(); ++j)
			{
				std::cout << Answ_coord[i][j] << ", ";
			}
			std::cout << std::endl;
		}


		auto zigota = pMandC;
		for(size_t i = 0; i < table.size(); ++i)
		{
			for(size_t j = 0; j < table.size(); ++j)
			{
				//scalar multuplication! this is a mistake
				//for(size_t beg = 0; beg < table.size(); ++beg)
				//{
				//	zigota[i][j] += (L_matrix[i][beg] * pMandC[beg][i]); 
				//}
				zigota[i][j] = L_matrix[i][j] * pMandC[j][i];
				std::cout << L_matrix[i][j] << ", ";
			}
			std::cout << std::endl;
		}
	
		for(size_t i = 0; i < table.size(); ++i)
		{
			for(size_t j = 0; j < table.size(); ++j)
			{
				std::cout << zigota[i][j] << ", ";
			}
			std::cout << std::endl;
		}


		float sum = 0.0;
		for(auto it = zigota.begin(); it != zigota.end(); ++it)
		{
			sum += std::accumulate((*it).begin(), (*it).end(), 0.0);			
			std::cout <<"sum " <<   sum << std::endl;
		}
		std::cout << "VICTORYYYYYYYYYYYYY : " << sum << std::endl;



		std::cout << "VIT_S " << vitrati << std::endl;
		std::cout << "1.0 - VIT_S " << 1.0 -   vitrati << std::endl;



	return ret;
}

int main()
{

	std::vector<float> pM(0, 0);
	std::vector<float> pK(0, 0);
	std::vector<std::vector<uint32_t>> table(20, std::vector<uint32_t>(20,0));
	
	std::vector<float> pC(20,0);
	std::vector<std::vector<float>> pMandC(20, std::vector<float>(20,0));
	std::vector<std::vector<float>> pMifC(20, std::vector<float>(20,0));
	


	collect_the_csv_values(pM, pK, pC, table);
	create_pMifC_csv(pM, pK, pC, table, pMandC, pMifC);
	auto answ_D = Deterministic(table, pMandC, pMifC);
	std::cout << std::endl << "Deterministic Function" << std::endl;
	for(auto it = answ_D.begin(); it != answ_D.end(); ++it)
	{
		std::cout << (*it).first << " ----> " << (*it).second << std::endl;
	}

	std::cout << std::endl << std::endl;
	auto answ_S = Stohastic(table, pMandC, pMifC);
	std::cout << std::endl << "Stohastic Function" << std::endl;
	for(auto it = answ_S.begin(); it != answ_S.end(); ++it)
	{
		std::cout << (*it).first << " ----> " << (*it).second << std::endl;
	}



return 0;
}
