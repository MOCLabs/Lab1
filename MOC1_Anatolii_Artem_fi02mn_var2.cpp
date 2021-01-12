#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <fstream>


void collect_the_csv_values(std::vector<float>& pM, std::vector<float>& pK, std::vector<float>& pC, std::vector<std::vector<uint32_t>>& table)
{
	std::fstream csv_read_probs("/home/ftl/moc1/vars/prob_02.csv", std::ios::in);
	std::fstream csv_read_table("/home/ftl/moc1/vars/table_02.csv", std::ios::in);

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
		}
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
	return CT_and_PT;
}

std::vector<std::pair<uint32_t, uint32_t>> Stohastic(std::vector<std::vector<uint32_t>> table, std::vector<std::vector<float>> pMandC, std::vector<std::vector<float>> pMifC)
{
	std::vector<std::pair<uint32_t, uint32_t>> ret(20, std::pair<uint32_t, uint32_t>(0,0));
	float vitrati = 0.0;
	uint32_t answ = 0;
	float L = 0.0;
	for(size_t col = 0; col < table.size(); ++col)
	{
		ret[col].first = col;
		answ = 0;
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

		L = 0.0;
		if(table[col][it] != it)
		{
			L = ((1.0) / (float)(amount));
		}
		vitrati += L * pMifC[it][col];
		ret[col].second = it;



		
	std::cout << "prob for " << col << " is " << ((float)(1.0) / ((float)amount)) << std::endl;
	std::cout << "vitrati for " << col << " is " << vitrati << std::endl;


	}	
	



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
