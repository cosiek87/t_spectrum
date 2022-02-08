////////////////////////////////////////////////////////////////////////////
//
//	Program do analizy plikow zawierajacych cale sekwencje pomiarow przeprowadzonych przy uzyciu ukladu
//  Caen DT5730 oraz programu Compass. Skrypt ten w glownej mierze sluzy tylko do pomiarow z 09.2021
//  Jednakze rdzen sortowania eventow i dopasowywania funckji jest ten dla wszystkich pomiarow
//  w ktorych uzywany byl uklad firmy Caen
//
//	Autor: Przemyslaw Sekowski
//
//	email: przemyslaw.sekowski@fuw.edu.pl
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <algorithm> 


UShort_t energia; // zmienne do ktorych wpadaja wartosci z kolejnych eventow
UShort_t energia_other;
ULong64_t czas;
UShort_t channel;
vector<ULong64_t> czas_koincydencji(3);
czas_koincydencji[0] = 100000, czas_koincydencji[1] = 16000, czas_koincydencji[2] = 40000;
ULong64_t closest_time_high, closest_time_low, delta_time_high, delta_time_low, delta_time;
Int_t nentries;						   // calkowita liczba eventow w pliku root
Int_t n_entried_entries = 0;		   //liczba eventow ktora przeszla analize/sortowanie
//wektory zawierajace zliczenia w pikach, ich bledy, oraz wartosci odcinania widma
vector<vector<Double_t>> zliczenia(6, vector<Double_t>(liczba_pomiarow));
vector<vector<Double_t>> blad_zliczenia(6, vector<Double_t>(liczba_pomiarow));
vector<vector<Double_t>> blad_zliczenia_fixed(6, vector<Double_t>(liczba_pomiarow));
vector<vector<Double_t>> zakres_energii(2, vector<Double_t>(6));
vector<vector<Double_t>> eliptyczny_zakres(2, vector<Double_t>(6));
TF1 *pre_dopasowanie[6]; //tworzona jest macierz funkcji dla kazdego widma z kazdej pozycji i detektora	
vector<ULong64_t> czasy;
vector<vector<vector<ULong64_t>>> wektor_timestamp(3, vector<vector<ULong64_t>>(liczba_pomiarow, czasy));
vector<Int_t> pozycje;
vector<vector<vector<Int_t>>> wektor_entry(3, vector<vector<Int_t>>(liczba_pomiarow, pozycje));
//zmienne potrzebne do okreslenia danego stanu tarczy (w trakcie obracania czy nie) i pozycje
Int_t stan, pozycja;
//wektor do ktore zapisane beda momenty rozpoczecia i zakonczenia pomiaru
vector<ULong64_t> wektor_czasu;
int delta_ch;
zakres_energii[0][0] = 900, zakres_energii[0][1] = 1800, zakres_energii[0][2] = 2400, zakres_energii[0][3] = 1650, zakres_energii[0][4] = 2050, zakres_energii[0][5] = 2400;  //wartosci energii wyrazonej w kanalach do dopasowania f. gaussa do widma
zakres_energii[1][0] = 1450, zakres_energii[1][1] = 2400, zakres_energii[1][2] = 2850, zakres_energii[1][3] = 1950, zakres_energii[1][4] = 2500, zakres_energii[1][5] = 3200; //wartosci energii wyrazonej w kanalach do dopasowania f. gaussa do widma
auto h_time = new TH1F("spek_time", "Widmo czasowe", 1e5, 0, 6e5);
auto h_step_time = new TH1F("spek_step_time", "Widmo krokow", 1e5, 0, 6e5);
auto h_rot_time = new TH1F("spek_rot_time", "Widmo momentow obrotu", 1e5, 0, 6e5);
auto h_delta_time = new TH1F("spek_delta_time", "Widmo delta time", 5e3, 0, 1e7);




//Funkcja będąca sumą funkcji Gaussa oraz f. liniowej. Paramtery: par[0] - liczba zliczen, par[1] - centroida, par[2] - sigma, par[3] - współczynnik kierunkowy, par[4] - wyraz wolny;
Double_t gausswithlinearbkg(Double_t *xarg, Double_t *par)
{
	Double_t x = xarg[0], result = 0.;
	result = par[0] / (par[2] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp((-TMath::Power((x - par[1]), 2)) / (2 * TMath::Power(par[2], 2))) + par[3] * x + par[4];
	return result;
}



//Funkcji Gaussa. Paramtery: par[0] - liczba zliczen, par[1] - centroida, par[2] - sigma;
Double_t mgauss(Double_t *xarg, Double_t *par)
{
	Double_t x = xarg[0], result = 0.;
	result = par[0] / (par[2] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp((-0.5 * TMath::Power((x - par[1]), 2)) / (TMath::Power(par[2], 2)));
	return result;
}



//Przybliżenie aglebraiczne wydajnosci pomiaru, opisane w ,,Evaluation of the Influence of Neighboring Radioactive Sources Placed on a Rotating Disk on the Photon Energy Spectrum'', T. Matulewicz et al.
Double_t przyblizenie_Mat(Double_t stosunek_r_d, Double_t pozycja, Double_t wszystkich_pozycji)
{
	Double_t result;
	result = 1 / (2 * pow(stosunek_r_d, 2) * (1 - cos(2.0 * TMath::Pi() * pozycja / wszystkich_pozycji)) + 1);
	return result;
}



//Funkcja zwracajaca numer pomiaru do ktorego nalezy event wzgledem podanego wektora czasu oraz czasu eventu;
UShort_t numer_pomiaru(vector<ULong64_t> vect, ULong64_t czas)
{
	UShort_t pozycja;
	for (Int_t i = 0; i < vect.size(); i++){
		pozycja = i;
		if (czas < vect[i]) break;
	}
	return pozycja;
}



ULong64_t search_closest(std::vector<ULong64_t>& sorted_array, ULong64_t x) {

    auto iter_geq = std::lower_bound(sorted_array.begin(), sorted_array.end(), x);
    if (iter_geq == sorted_array.begin()) return 0;
    ULong64_t a = *(iter_geq - 1);
    ULong64_t b = *(iter_geq);
    if (fabs(x - a) < fabs(x - b)) return iter_geq - sorted_array.begin() - 1;
    return iter_geq - sorted_array.begin();
}



//Funkcja zwracajaca czy energia z danego eventu zawiera sie przedziale energetycznym 511 keV +/- 20%
//Nie jest to funkcja obowiazkowa
Bool_t w_zakresie_511kev(UShort_t energia, UShort_t channel)
{
	Bool_t w_zakresie = true;
	switch (channel){
	case 0:
		if (energia > 800 && energia < 1500) w_zakresie = false;
		break;
	case 1:
		if (energia > 1700 && energia < 2500) w_zakresie = false;
		break;
	case 2:
		if (energia > 2200 && energia < 3000) w_zakresie = false;
		break;
	case 3:
		if (energia > 1400 && energia < 2200) w_zakresie = false;
		break;
	case 4:
		if (energia > 1800 && energia < 2600) w_zakresie = false;
		break;
	case 5:
		if (energia > 2400 && energia < 3200) w_zakresie = false;
		break;
	}
	return w_zakresie;
}



//Funkcja oblicza wektor czasu na podstawie zbioru eventow pochodzacych od silnika. Estymacja momentow obrotu polega na klasteryzacji krokow a nastepnie wybieranie skrajnych czasow z klastra;
//Do funkcji tej nalezy wpisac nazwe pliku z danymi (SDataR_*.root)
//Zwraca wektor czasu, ktory sluzy do oszacowania numeru pomiaru
std::vector< ULong64_t > obliczanie_wektora_czasu(std::string fileName)
{
 	std::ifstream in(fileName.c_str());
 	std::vector<ULong64_t> vecOfStr; 		// Sprawdza czy plik jest ok
	std::string str; 						// Czyta kolejne linijki dopoki plik sie nie skonczy
    while (std::getline(in, str)){ 			// Jesli string jest niezerowej dlugosci to jest zapisywany do wektora czasu
        if(str.size() > 0) vecOfStr.push_back(stod(str)*1e12);
    }
    in.close(); 							//Zamyka plik
	return vecOfStr;						// Zwraca wektor czasu
}



void wektory_czasu_koincydencji(int liczba_pomiarow, std::vector<ULong64_t> wektor_czasu,std::vector<std::vector<std::vector<ULong64_t>>>& wektor_timestamp, std::vector<std::vector<std::vector<Int_t>>>& wektor_entry)
{
	cout<<"tworzy sie macierz koincydencji"<<endl;
	Char_t nazwa_1[200] = "SDataR_BN_100Gy_1.root";
	ULong64_t czas;
	UShort_t channel, stan, pozycja, energia;
	Int_t n_times;
	auto *f_1 = new TFile(nazwa_1);
	auto *t_1 = (TTree *)f_1->Get("Data_R");
	t_1->SetBranchAddress("Timestamp", &czas);
	t_1->SetBranchAddress("Channel", &channel);
	t_1->SetBranchAddress("Energy", &energia);
	n_times = (Int_t)t_1->GetEntries();
	for (Int_t i = 0; i < n_times; i++){
		t_1->GetEntry(i);
		if (w_zakresie_511kev(energia, channel)) continue;
		stan = numer_pomiaru(wektor_czasu, czas); 						//stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
		pozycja = stan / 2;
		if (channel > 5 && stan % 2 == 0) continue;
		if ( channel % 2 == 0) continue;
		wektor_timestamp[(channel-1)/2][pozycja].push_back(czas); 		// wypelniany jest wektor z czasem do koincydencji
		wektor_entry[(channel-1)/2][pozycja].push_back(i); 				// wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
																		// if (i%10000==0) cout<<"wrzucono do "<<channel/2<<" oraz pozycji "<<pozycja<<" wartosci "<<czas<<" oraz "<<i<<endl;
	}
	cout<<"stworzyla sie macierz koincydencji"<<endl;
}



ULong64_t binary_search(std::vector<ULong64_t> wektor, ULong64_t value)
{
return (TMath::BinarySearch(std::begin(wektor), std::end(wektor), value) - wektor.begin());
}



ULong64_t closest(std::vector<ULong64_t> const& vec, ULong64_t value) {
    auto const it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end()) { return vec.size() - 1; }
	if (it == vec.begin()) { return 0; }
    return (it-vec.begin());
}


void w_zakresie_elipsy(UShort_t energia_dol, UShort_t energia_gora, UShort_t channel, TF1 *pre_dopasowanie, Double_t ile_sigma) {
	UShort_t srodek_1 = pre_dopasowanie[int(channel/2)]->GetParameter(1), srodek_2 = pre_dopasowanie[int(channel/2)+1]->GetParameter(1);
	UShort_t sigma_1 = pre_dopasowanie[int(channel/2)]->GetParameter(2), sigma_2 = pre_dopasowanie[int(channel/2)+1]->GetParameter(2)
	if (((energia_dol-srodek_1)/sigma_1)**2 + ((energia_gora-srodek_2)/sigma_2)**2 > 1) return false;
	return true;
}



void creating_histograms(TH1F *h, TH1F *total_h, TH2F *h2d, TH2F *total_h2d)
{
	char name[20];
	char title[100];
	for (Int_t i = 0; i < 6; i++)	{
		for (Int_t m = 0; m < liczba_pomiarow; m++){
			sprintf(name, "spek_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum pozycji h%d det %d", m + 1, i);
			h[i][m] = new TH1F(name, title, numberofbins, minimum, maksimum);
			sprintf(name, "spek_2D_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum 2D pozycji h%d det %d", m + 1, i);
			if (i%2 == 0) h_2d[i/2][m] = new TH2F(name, title, 150, 500, 3500, 150, 500, 3500);
		}
		sprintf(name, "spek_total_det_%d", i);
		sprintf(title, "Spektrum calkowite det %d", i);
		total_h[i] = new TH1F(name, title, numberofbins*4, minimum, maksimum);
	}
	for (Int_t i = 0; i < 3; i++){
		sprintf(name, "spek_total_2d_det_%d", i);
		sprintf(title, "Spektrum calkowite 2D det %d", i);
		total_h_2d[i] = new TH2F(name, title, numberofbins*4, minimum, maksimum, numberofbins*4, minimum, maksimum);
	}
	cout << "Zakonczono tworzenie histogramow" <<  endl;
}



void filling_total_histograms(TTree *t_1, Int_t nentries, std::vector<ULong64_t> wektor_czasu, TH1F *h_step_time, TH1F *total_h, Int_t para_detektorow)
{
	for (Int_t i = 0; i < nentries; i++){
		t_1->GetEntry(i);
		if (w_zakresie_511kev(energia, channel)) continue;
		stan = numer_pomiaru(wektor_czasu, czas); 				//stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
		if (channel == 7) h_step_time->Fill(czas / 1e12); 					// widmo krokow silnika
		pozycja = stan / 2;
		if (channel > 5 && stan % 2 == 0) continue;
		if (channel % 2 == 1) continue;	// jesli event jest zebrany na kanale innym niz do ktorych byly podlaczone detektory, kod przechodzi do kolejnego eventu
		total_h[channel]->Fill(energia);
	}
}


void filling_each_measurement_histograms(TTree *t_1, Int_t nentries, TH1F *h, TH2F *h_2d, std::vector<ULong64_t> wektor_czasu, TF1 *pre_dopasowanie, std::vector<ULong64_t> eliptic_par, Int_t para_detektorow){
	for (Int_t i = 0; i < nentries; i++){
		t_1->GetEntry(i);
		if (w_zakresie_511kev(energia, channel)) continue;
		stan = numer_pomiaru(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
		pozycja = stan / 2;
		if (channel > 5 && stan % 2 == 0) continue;	
		if (channel % 2 == 1) continue;	// jesli event jest zebrany na kanale innym niz do ktorych byly podlaczone detektory, kod przechodzi do kolejnego eventu
		auto closest_higher_time_index = closest(wektor_timestamp[(channel)/2][pozycja],czas);
		auto closest_lower_time_index = (closest_higher_time_index-1)*((closest_higher_time_index-1)>0);
		if (closest_higher_time_index > wektor_timestamp[(channel)/2][pozycja].size()-1 
		    || closest_lower_time_index > wektor_timestamp[(channel)/2][pozycja].size()-1 
			|| closest_higher_time_index > wektor_entry[(channel)/2][pozycja].size()-1 
			|| closest_lower_time_index > wektor_entry[(channel)/2][pozycja].size()-1) continue;
		closest_time_high = wektor_timestamp[(channel)/2][pozycja][closest_higher_time_index];
		closest_time_low = wektor_timestamp[(channel)/2][pozycja][closest_lower_time_index];
		delta_time_high = max(closest_time_high,czas) - min(closest_time_high,czas);
		delta_time_low = max(closest_time_low,czas) - min(closest_time_low,czas);
		auto closest_time_index = closest_higher_time_index;
		delta_time = delta_time_high;
		if (delta_time_high>delta_time_low){
			closest_time_index = closest_lower_time_index;
			delta_time = delta_time_low;
		}
		h_delta_time->Fill(delta_time);
		// cout<<wektor_timestamp[(channel-1)/2][pozycja][closest_time_index]<<" najblizszy czas do "<<czas<<endl;
		// cout << "kanal "<<channel<<endl;
		// cout << "delta time "<<delta_time<<endl;
		if (delta_time < czas_koincydencji[(channel)/2]){
			energia_other = energia;
			t_1->GetEntry(wektor_entry[(channel)/2][pozycja][closest_time_index]);
			if (w_zakresie_elipsy(energia_other, energia, channel, pre_dopasowanie, 3.5)) continue;
			h[channel][pozycja]->Fill(energia_other); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h[channel][pozycja]->Fill(energia); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h_2d[(channel)/2][pozycja]->Fill(energia_other,energia);
			total_h_2d[(channel)/2]->Fill(energia_other,energia);
			n_entried_entries++;				//zwiekszana jest liczba oznaczajaca eventy ktore przeszly analize
			}
		else continue;
	}
}



void fit_total_histograms(TF1 *pre_dopasowanie, TH1F *total_h, std::vector<Double_t> zakres_energii, Double_t zliczenia_amplituda, Int_t para_detektorow)
{
	for (Int_t det = 0; det < 6; det++){
			pre_dopasowanie[det] = new TF1("dopasowanie", gausswithlinearbkg, zakres_energii[0][det], zakres_energii[1][det], 5);
			pre_dopasowanie[det]->SetParameters(20000, (zakres_energii[0][det] + zakres_energii[1][det]) / 2, 25, -0.00001, 10);
			pre_dopasowanie[det]->SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
			pre_dopasowanie[det]->SetParLimits(0, 20*zliczenia_amplituda, 20000*zliczenia_amplituda);
			pre_dopasowanie[det]->SetParLimits(1, zakres_energii[0][det], zakres_energii[1][det]);
			pre_dopasowanie[det]->SetParLimits(2, 25, 100);
			pre_dopasowanie[det]->SetParLimits(3, -0.1, 0.1);
			pre_dopasowanie[det]->SetParLimits(4, -10000, 10000);
			total_h[det]->Fit("dopasowanie", "LMQR", "", zakres_energii[0][det], zakres_energii[1][det]);
			total_h[det]->Draw();
	}
}


void extract_number_of_counts(Int_t liczba_pomiarow, vector<vector<Double_t>> zliczenia, vector<vector<Double_t>> blad_zliczenia, Int_t para_detektorow)
{
	wektor_czasu.insert(wektor_czasu.begin(), 0); // na poczatek wektora czasu dorzucany jest moment rozpoczecia pomiaru.

	for (Int_t det = 0; det < 6; det++){
		cout << "Liczba zliczen w pomiarze - detektor " << det << endl;
		for (Int_t i = 0; i < liczba_pomiarow; i++) zliczenia[det][i] = (h[det][j]->IntegralAndError(0,17000,0,17000, blad_zliczenia[det][i], true,"width"));
	}

	for (Int_t det = 0; det < 6; det++){
		cout << "Liczba zliczen w pomiarze - detektor " << det << endl;
		for (Int_t i = 0; i < liczba_pomiarow; i++){
			cout << "\t\t" << zliczenia[det][i];
			cout << "\t\t" << blad_zliczenia[det][i]<< endl;
		}
	}

}