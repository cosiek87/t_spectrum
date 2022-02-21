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

using namespace std;

Int_t liczba_pomiarow = 42; // zmienna do okreslania wielkosci wektora w ktorym zapisywane sa zliczenia w piku, 
                            // odpowiada liczbie pojedynczych pomiarow
UShort_t energia; // zmienne do ktorych wpadaja wartosci z kolejnych eventow
UShort_t energia_other;
ULong64_t czas;
UShort_t channel;
vector < ULong64_t > czas_koincydencji(3);
ULong64_t closest_time_high, closest_time_low, delta_time_high, delta_time_low, delta_time;
Int_t nentries; // calkowita liczba eventow w pliku root
Int_t n_entried_entries = 0; //liczba eventow ktora przeszla analize/sortowanie
//wektory zawierajace zliczenia w pikach, ich bledy, oraz wartosci odcinania widma
vector < vector < Double_t >> zakres_energii(2, vector < Double_t > (6));
TF1 * pre_dopasowanie[6]; //tworzona jest macierz funkcji dla kazdego widma z kazdej pozycji i detektora	
vector < ULong64_t > czasy;

vector < Int_t > pozycje;

//zmienne potrzebne do okreslenia danego stanu tarczy (w trakcie obracania czy nie) i pozycje
Int_t stan, pozycja;
//wektor do ktore zapisane beda momenty rozpoczecia i zakonczenia pomiaru
vector < ULong64_t > wektor_czasu;
int delta_ch;
auto h_time = new TH1F("spek_time", "Widmo czasowe", 1e5, 0, 6e5);
auto h_step_time = new TH1F("spek_step_time", "Widmo krokow", 1e5, 0, 6e5);
auto h_rot_time = new TH1F("spek_rot_time", "Widmo momentow obrotu", 1e5, 0, 6e5);
auto h_delta_time = new TH1F("spek_delta_time", "Widmo delta time", 5e3, 0, 1e7);

// Funkcja będąca sumą funkcji Gaussa oraz f. liniowej. Paramtery: par[0] - liczba zliczen, 
// par[1] - centroida, par[2] - sigma, par[3] - współczynnik kierunkowy, par[4] - wyraz wolny;
Double_t gausswithlinearbkg(Double_t * xarg, Double_t * par) {
    Double_t x = xarg[0], result = 0.;
    result = par[0] / (par[2] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp((-TMath::Power((x - par[1]), 2)) 
            / (2 * TMath::Power(par[2], 2))) + par[3] * x + par[4];
    return result;
}

// Funkcji Gaussa. Paramtery: par[0] - liczba zliczen, par[1] - centroida, par[2] - sigma;
Double_t mgauss(Double_t * xarg, Double_t * par) {
    Double_t x = xarg[0], result = 0.;
    result = par[0] / (par[2] * TMath::Sqrt(2 * TMath::Pi())) * TMath::Exp((-0.5 * TMath::Power((x - par[1]), 2)) 
            / (TMath::Power(par[2], 2)));
    return result;
}

// Przybliżenie aglebraiczne wydajnosci pomiaru, opisane w 
// ,,Evaluation of the Influence of Neighboring Radioactive Sources Placed on a Rotating Disk on the Photon Energy Spectrum'', 
// T. Matulewicz et al.
Double_t przyblizenie_Mat(Double_t stosunek_r_d, Double_t pozycja, Double_t wszystkich_pozycji) {
    Double_t result;
    result = 1 / (2 * pow(stosunek_r_d, 2) * (1 - cos(2.0 * TMath::Pi() * pozycja / wszystkich_pozycji)) + 1);
    return result;
}

//Funkcja zwracajaca numer pomiaru do ktorego nalezy event wzgledem podanego wektora czasu oraz czasu eventu;
UShort_t numer_pomiaru(vector < ULong64_t > vect, ULong64_t czas) {
    UShort_t pozycja;
    for (Int_t i = 0; i < vect.size(); i++) {
        pozycja = i;
        if (czas < vect[i]) break;
    }
    return pozycja;
}

ULong64_t search_closest(std::vector < ULong64_t > & sorted_array, ULong64_t x) {

    auto iter_geq = std::lower_bound(sorted_array.begin(), sorted_array.end(), x);
    if (iter_geq == sorted_array.begin()) return 0;
    ULong64_t a = * (iter_geq - 1);
    ULong64_t b = * (iter_geq);
    if (fabs(x - a) < fabs(x - b)) return iter_geq - sorted_array.begin() - 1;
    return iter_geq - sorted_array.begin();
}

//Funkcja zwracajaca czy energia z danego eventu zawiera sie przedziale energetycznym 511 keV +/- 20%
//Nie jest to funkcja obowiazkowa
Bool_t w_zakresie_511kev(UShort_t energia, UShort_t channel) {
    Bool_t w_zakresie = true;
    switch (channel) {
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

// Funkcja oblicza wektor czasu na podstawie zbioru eventow pochodzacych od silnika. 
// Estymacja momentow obrotu polega na klasteryzacji krokow a nastepnie wybieranie skrajnych czasow z klastra;
// Do funkcji tej nalezy wpisac nazwe pliku z danymi (SDataR_*.root)
// Zwraca wektor czasu, ktory sluzy do oszacowania numeru pomiaru
std::vector < ULong64_t > obliczanie_wektora_czasu(std::string fileName) {
    std::ifstream in (fileName.c_str());
    std::vector < ULong64_t > vecOfStr; // Sprawdza czy plik jest ok
    std::string str; // Czyta kolejne linijki dopoki plik sie nie skonczy
    while (std::getline( in , str)) { // Jesli string jest niezerowej dlugosci to jest zapisywany do wektora czasu
        if (str.size() > 0) vecOfStr.push_back(stod(str) * 1e12);
    } in .close(); //Zamyka plik
    return vecOfStr; // Zwraca wektor czasu
}

void wektory_czasu_koincydencji(int liczba_pomiarow, std::vector < ULong64_t > wektor_czasu,
                                std::vector < std::vector < std::vector < ULong64_t >>> & wektor_timestamp,
                                std::vector < std::vector < std::vector < Int_t >>> & wektor_entry,
                                int para_detektorow) {
    Int_t liczba_par_det, start_det, stop_det, liczba_det;
	if (para_detektorow == -1) {
		liczba_par_det = 3;
		liczba_det = 6;
		start_det = 0;
		stop_det = 6;
	}
	else {
		liczba_par_det = 1;
		liczba_det = 2;
		start_det = 2 * para_detektorow;
		stop_det = start_det + 2;
	}
    Char_t nazwa_1[200] = "E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\SDataR_BN_100Gy_1.root";
    ULong64_t czas;
    UShort_t channel, stan, pozycja, energia;
    Int_t n_times;
    auto * f_1 = new TFile(nazwa_1);
    auto * t_1 = (TTree * ) f_1 -> Get("Data_R");
    t_1 -> SetBranchAddress("Timestamp", & czas);
    t_1 -> SetBranchAddress("Channel", & channel);
    t_1 -> SetBranchAddress("Energy", & energia);
    n_times = (Int_t) t_1 -> GetEntries();
    for (Int_t i = 0; i < n_times; i++) {
        t_1 -> GetEntry(i);
        if (w_zakresie_511kev(energia, channel)) continue;
        stan = numer_pomiaru(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
        pozycja = stan / 2;
        if (channel > 5 && stan % 2 == 0) continue;
        if (channel % 2 == 0) continue;
    	if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
        wektor_timestamp[(channel - 1) / 2][pozycja].push_back(czas); // wypelniany jest wektor z czasem do koincydencji
        wektor_entry[(channel - 1) / 2][pozycja].push_back(i); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
        // if (i%10000==0) cout<<"wrzucono do "<<channel/2<<" oraz pozycji "<<pozycja<<" wartosci "<<czas<<" oraz "<<i<<endl;
    }
    cout << "stworzyla sie macierz koincydencji" << endl;
}

ULong64_t binary_search(std::vector < ULong64_t > wektor, ULong64_t value) {
    return (TMath::BinarySearch(std::begin(wektor), std::end(wektor), value) - wektor.begin());
}

ULong64_t closest(std::vector < ULong64_t >
    const & vec, ULong64_t value) {
    auto
    const it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end()) {
        return vec.size() - 1;
    }
    if (it == vec.begin()) {
        return 0;
    }
    return (it - vec.begin());
}

bool w_zakresie_elipsy(UShort_t energia_dol, UShort_t energia_gora, UShort_t channel,TF1 * pre_dopasowanie[], Double_t ile_sigma) {
    UShort_t srodek_1 = pre_dopasowanie[channel-1] -> GetParameter(1), srodek_2 = pre_dopasowanie[channel] -> GetParameter(1);
    UShort_t sigma_1 = pre_dopasowanie[channel-1] -> GetParameter(2), sigma_2 = pre_dopasowanie[channel] -> GetParameter(2);
    if (TMath::Power((energia_dol - srodek_1) / (ile_sigma*sigma_1), 2) + TMath::Power((energia_gora - srodek_2) 
        / (ile_sigma*sigma_2), 2) > 1) return true;
    return false;
}
