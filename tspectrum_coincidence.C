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
    Int_t para_detektorow = 2;
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
    Int_t numberofbins = 1700; //liczba binow w widme energergetycznym; zmniejszajac jej wartosc otrzymujemy lepsza statystykew binach ale gorsza rozdzielczosc (mniej pewna wartosc dopasowanej centroidy)
    Float_t minimum = 0, maksimum = 17000; //zakres widma energetycznego wyrazone w kanalach
    Double_t zliczenia_amplituda = maksimum / numberofbins;
    TH1F * h[liczba_det][liczba_pomiarow];
    TH1F * total_h[liczba_det];
    TH2F * h_2d[liczba_par_det][liczba_pomiarow];
    TH2F * total_h_2d[liczba_par_det];
    TF1 * dopasowanie[liczba_det][liczba_pomiarow]; //tworzona jest macierz funkcji dla kazdego widma z kazdej pozycji i detektora
    Char_t nazwa_1[200] = "E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\SDataR_BN_100Gy_1.root"; //plik zrodlowy z danymi
    auto * f_1 = new TFile(nazwa_1); //tworzone sa tu plik wyjsciowy i wszystkie histogramy, a takze do zmiennych energia, czas i channel przypisywane sa wartosci eventow
    auto * t_1 = (TTree * ) f_1 -> Get("Data_R");
    t_1 -> SetBranchAddress("Energy", & energia);
    t_1 -> SetBranchAddress("Timestamp", & czas);
    t_1 -> SetBranchAddress("Channel", & channel);
    auto * f_output = new TFile("E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\Analysis_SDataR_BN_100Gy_1_koincydencja_w_total.root", "RECREATE");
    vector < vector < Double_t >> zliczenia(liczba_det, vector < Double_t > (liczba_pomiarow));
    vector < vector < Double_t >> blad_zliczenia(liczba_det, vector < Double_t > (liczba_pomiarow));
    vector < vector < Double_t >> blad_zliczenia_fixed(liczba_det, vector < Double_t > (liczba_pomiarow));
    vector < vector < vector < ULong64_t >>> wektor_timestamp(liczba_par_det, vector < vector < ULong64_t >> (liczba_pomiarow, czasy));
    vector < vector < vector < Int_t >>> wektor_entry(liczba_par_det, vector < vector < Int_t >> (liczba_pomiarow, pozycje));

    zakres_energii[0][0] = 900, zakres_energii[0][1] = 1800, zakres_energii[0][2] = 2400, zakres_energii[0][3] = 1650, zakres_energii[0][4] = 2050, zakres_energii[0][5] = 2400; //wartosci energii wyrazonej w kanalach do dopasowania f. gaussa do widma
    zakres_energii[1][0] = 1450, zakres_energii[1][1] = 2400, zakres_energii[1][2] = 2850, zakres_energii[1][3] = 1950, zakres_energii[1][4] = 2500, zakres_energii[1][5] = 3200; //wartosci energii wyrazonej w kanalach do dopasowania f. gaussa do widma
	cout<<"Tworzenie histogramow"<<endl;
    char name[20];
	char title[100];


	for (Int_t i = 0; i < liczba_det; i++) {
		for (Int_t m = 0; m < liczba_pomiarow; m++) {
			sprintf(name, "spek_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum pozycji h%d det %d", m + 1, i);
			h[i][m] = new TH1F(name, title, numberofbins, minimum, maksimum);
			sprintf(name, "spek_2D_%d_det_%d", m + 1, i);
			sprintf(title, "Spektrum 2D pozycji h%d det %d", m + 1, i);
			if (i % 2 == 0) h_2d[i / 2][m] = new TH2F(name, title, 150, 500, 3500, 150, 500, 3500);
		}
		cout<<"total det: "<<i<<endl;
		sprintf(name, "spek_total_det_%i", i);
		sprintf(title, "Spektrum calkowite det %i", i);
		total_h[i] = new TH1F(name, title, 4000, minimum, maksimum);
		cout<<"Utworzono histogram: name: "<<name<<" title: "<<title<<" obiekt "<<&total_h[i]<<endl;
	}
	cout<<"Przygotowywanie widm czastkowych 2d"<<endl;
	for (Int_t i = 0; i < liczba_par_det; i++) {
		sprintf(name, "spek_total_2d_det_%d", i);
		sprintf(title, "Spektrum calkowite 2D det %d", i);
		total_h_2d[i] = new TH2F(name, title, 400, minimum, maksimum, 400, minimum, maksimum);
	}
	cout << "Zakonczono tworzenie histogramow" << endl;

	cout<<"Obliczanie wektora czasu"<<endl;
    wektor_czasu = obliczanie_wektora_czasu("E:\\EKSPERYMENT\\09.2021\\21.09\\BN\\DAQ\\BN_100Gy_1\\RAW\\os_czasu_bn.txt");
    cout << "wektor czasu: " << endl;
    for (Int_t k = 0; k < wektor_czasu.size(); k++) {
        h_rot_time -> Fill(wektor_czasu[k] / 1e12); //wypelniany jest histogram momentow obrotu tarcza
    }
    nentries = (Int_t) t_1 -> GetEntries();


	cout<<"Tworzenie wektorow koincydencji"<<endl;
    wektory_czasu_koincydencji(liczba_pomiarow, wektor_czasu, wektor_timestamp, wektor_entry, para_detektorow);


	cout<<"Wypelnianie histogramow"<<endl;
    for (Int_t i = 0; i < nentries; i++) {
		t_1 -> GetEntry(i);
		if (w_zakresie_511kev(energia, channel)) continue;
		stan = numer_pomiaru(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
		if (channel == 7) h_step_time -> Fill(czas / 1e12); // widmo krokow silnika
		pozycja = stan / 2;
		if (channel > 5 && stan % 2 == 0) continue;
		if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
		total_h[channel] -> Fill(energia);
	}
	
	for (Int_t det = 0; det < liczba_det; det++) {
		pre_dopasowanie[det] = new TF1("dopasowanie", gausswithlinearbkg, zakres_energii[0][start_det+det], zakres_energii[1][start_det+det], 5);
		pre_dopasowanie[det] -> SetParameters(20000, (zakres_energii[0][start_det+det] + zakres_energii[1][start_det+det]) / 2, 25, -0.00001, 10);
		pre_dopasowanie[det] -> SetParNames("Amplituda", "Srednia", "Sigma", "A", "B");
		pre_dopasowanie[det] -> SetParLimits(0, 20 * zliczenia_amplituda, 200000 * zliczenia_amplituda);
		pre_dopasowanie[det] -> SetParLimits(1, zakres_energii[0][start_det+det], zakres_energii[1][start_det+det]);
		pre_dopasowanie[det] -> SetParLimits(2, 25, 100);
		pre_dopasowanie[det] -> SetParLimits(3, -0.1, 0.1);
		pre_dopasowanie[det] -> SetParLimits(4, -10000, 10000);
		total_h[det] -> Fit("dopasowanie", "LMQR0", "", zakres_energii[0][start_det+det], zakres_energii[1][start_det+det]);
		// total_h[det] -> Draw();

	}
	for (int i = 0 ; i<liczba_det; i++) {
		cout<<" obiekt "<<&total_h[i]<<endl;
		cout<<total_h[i] -> GetEntries()<<endl;
		cout<<zakres_energii[0][i]<<" " <<zakres_energii[1][i]<<endl;
		cout<< "centroida "<<pre_dopasowanie[i] -> GetParameter(1)<<" sigma "<< pre_dopasowanie[i] -> GetParameter(2)<<endl;;
	}
	cout<<"Wypelnianie widm kazdego pomiaru"<<endl;
    for (Int_t i = 0; i < nentries; i++) {
		t_1 -> GetEntry(i);
		if (w_zakresie_511kev(energia, channel)) continue;
		stan = numer_pomiaru(wektor_czasu, czas); //stan tarczy: wartosc nieparzysta - obrot, parzysta - pomiar
		pozycja = stan / 2;
		if (channel > 5 && stan % 2 == 0) continue;
		if (channel % 2 == 1) continue; // jesli event jest zebrany na kanale innym niz do ktorych byly podlaczone detektory, kod przechodzi do kolejnego eventu
		if (channel < start_det || channel >= stop_det) continue;
		channel %= liczba_det;
		auto closest_higher_time_index = closest(wektor_timestamp[(channel) / 2][pozycja], czas);
		auto closest_lower_time_index = (closest_higher_time_index - 1) * ((closest_higher_time_index - 1) > 0);
		if (closest_higher_time_index > wektor_timestamp[(channel) / 2][pozycja].size() - 1 ||
			closest_lower_time_index > wektor_timestamp[(channel) / 2][pozycja].size() - 1 ||
			closest_higher_time_index > wektor_entry[(channel) / 2][pozycja].size() - 1 ||
			closest_lower_time_index > wektor_entry[(channel) / 2][pozycja].size() - 1) continue;
		closest_time_high = wektor_timestamp[(channel) / 2][pozycja][closest_higher_time_index];
		closest_time_low = wektor_timestamp[(channel) / 2][pozycja][closest_lower_time_index];
		delta_time_high = max(closest_time_high, czas) - min(closest_time_high, czas);
		delta_time_low = max(closest_time_low, czas) - min(closest_time_low, czas);
		auto closest_time_index = closest_higher_time_index;
		delta_time = delta_time_high;
		if (delta_time_high > delta_time_low) {
			closest_time_index = closest_lower_time_index;
			delta_time = delta_time_low;
		}	
		h_delta_time -> Fill(delta_time);
		// cout<<wektor_timestamp[(channel-1)/2][pozycja][closest_time_index]<<" najblizszy czas do "<<czas<<endl;
		// cout << "kanal "<<channel<<endl;
		// cout << "delta time "<<delta_time<<endl;
		if (i%10000 == 0) cout<<i<<endl;
		if (delta_time < czas_koincydencji[(channel) / 2]) {
			energia_other = energia;
			t_1 -> GetEntry(wektor_entry[(channel) / 2][pozycja][closest_time_index]);
			if (channel < start_det || channel >= stop_det) continue;
			channel %= liczba_det;	
			if (w_zakresie_elipsy(energia_other, energia, channel, pre_dopasowanie, 3)) continue;
			h[channel-1][pozycja] -> Fill(energia_other); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h[channel][pozycja] -> Fill(energia); // wypelniane jest widmo energetyczne w zaleznosci od kanalu i pozycji
			h_2d[(channel) / 2][pozycja] -> Fill(energia_other, energia);
			total_h_2d[(channel) / 2] -> Fill(energia_other, energia);
			n_entried_entries++; //zwiekszana jest liczba oznaczajaca eventy ktore przeszly analize
		} else continue;
	}
	cout<<"Eksport danych"<<endl;
    wektor_czasu.insert(wektor_czasu.begin(), 0); // na poczatek wektora czasu dorzucany jest moment rozpoczecia pomiaru.

	for (Int_t det = 0; det < liczba_det; det++) {
		for (Int_t i = 0; i < liczba_pomiarow; i++) zliczenia[det][i] = (h[det][i]->IntegralAndError(0, 17000, blad_zliczenia[det][i]));
	}

	for (Int_t det = 0; det < liczba_det; det++) {
		cout << "Liczba zliczen w pomiarze - detektor " << start_det + det << endl;
		for (Int_t i = 0; i < liczba_pomiarow; i++) {
			cout << "\t\t" << zliczenia[det][i];
			cout << "\t\t" << blad_zliczenia[det][i] << endl;
		}
	}

    f_output -> Write();
}