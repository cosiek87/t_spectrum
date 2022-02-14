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

#include "tspectrum_coincidence.hpp"

using namespace std;

//Glowna czesc skryptu
void tspectrum_coincidence() {

    czas_koincydencji = {
        100000,
        16000,
        40000
    };
    Int_t liczba_pomiarow = 42; // zmienna do okreslania wielkosci wektora w ktorym zapisywane sa zliczenia w piku, odpowiada liczbie pojedynczych pomiarow
    Int_t para_detektorow = -1;
    Int_t numberofbins = 425; //liczba binow w widme energergetycznym; zmniejszajac jej wartosc otrzymujemy lepsza statystykew binach ale gorsza rozdzielczosc (mniej pewna wartosc dopasowanej centroidy)
    Float_t minimum = 0, maksimum = 17000; //zakres widma energetycznego wyrazone w kanalach
    Double_t zliczenia_amplituda = maksimum / numberofbins;
    TH1F * h[6][liczba_pomiarow];
    TH1F * total_h[6];
    TH2F * h_2d[3][liczba_pomiarow];
    TH2F * total_h_2d[3];
    TF1 * dopasowanie[6][liczba_pomiarow]; //tworzona jest macierz funkcji dla kazdego widma z kazdej pozycji i detektora
    Char_t nazwa_1[200] = "E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\SDataR_BN_100Gy_1.root"; //plik zrodlowy z danymi
    auto * f_1 = new TFile(nazwa_1); //tworzone sa tu plik wyjsciowy i wszystkie histogramy, a takze do zmiennych energia, czas i channel przypisywane sa wartosci eventow
    auto * t_1 = (TTree * ) f_1 -> Get("Data_R");
    t_1 -> SetBranchAddress("Energy", & energia);
    t_1 -> SetBranchAddress("Timestamp", & czas);
    t_1 -> SetBranchAddress("Channel", & channel);
    auto * f_output = new TFile("E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\Analysis_SDataR_BN_100Gy_1_koincydencja_w_total.root", "RECREATE");
    vector < vector < Double_t >> zliczenia(6, vector < Double_t > (liczba_pomiarow));
    vector < vector < Double_t >> blad_zliczenia(6, vector < Double_t > (liczba_pomiarow));
    vector < vector < Double_t >> blad_zliczenia_fixed(6, vector < Double_t > (liczba_pomiarow));
    vector < vector < vector < ULong64_t >>> wektor_timestamp(3, vector < vector < ULong64_t >> (liczba_pomiarow, czasy));
    vector < vector < vector < Int_t >>> wektor_entry(3, vector < vector < Int_t >> (liczba_pomiarow, pozycje));
    zakres_energii[0][0] = 900, zakres_energii[0][1] = 1800, zakres_energii[0][2] = 2400, zakres_energii[0][3] = 1650, zakres_energii[0][4] = 2050, zakres_energii[0][5] = 2400; //wartosci energii wyrazonej w kanalach do dopasowania f. gaussa do widma
    zakres_energii[1][0] = 1450, zakres_energii[1][1] = 2400, zakres_energii[1][2] = 2850, zakres_energii[1][3] = 1950, zakres_energii[1][4] = 2500, zakres_energii[1][5] = 3200; //wartosci energii wyrazonej w kanalach do dopasowania f. gaussa do widma

    creating_histograms(h, total_h, h_2d, total_h_2d, para_detektorow);

    wektor_czasu = obliczanie_wektora_czasu("E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\os_czasu_bn.txt");
    cout << "wektor czasu: " << endl;
    for (Int_t k = 0; k < wektor_czasu.size(); k++) {
        h_rot_time -> Fill(wektor_czasu[k] / 1e12); //wypelniany jest histogram momentow obrotu tarcza
    }
    nentries = (Int_t) t_1 -> GetEntries();
    wektory_czasu_koincydencji(liczba_pomiarow, wektor_czasu, wektor_timestamp, wektor_entry, para_detektorow);
    filling_total_histograms(t_1, wektor_czasu, h_step_time, total_h, para_detektorow);
    f_output -> Write(); //zapisywane sa widma w pliku wyjsciowym
    fit_total_histograms(pre_dopasowanie, total_h, zakres_energii, zliczenia_amplituda, para_detektorow); //do kazdego widma dopasowana zostaje f. guassa z liniowym tlem
    filling_each_measurement_histograms(t_1, wektor_czasu, h_step_time, total_h, pre_dopasowanie, para_detektorow);
    extract_number_of_counts();

    f_output -> Write();
}