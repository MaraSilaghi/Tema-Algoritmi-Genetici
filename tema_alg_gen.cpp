//Input fisier:
// -1 1 2
// -1 2
// 6
// 20
// 0.25
// 0.01
// 50

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <bitset>
#include <iomanip>
#include <fstream>
#include <climits>

using namespace std;
ifstream fin("input.in");
ofstream fout("evolutie.txt");

double a, b, c, limInf, limSup, probCrossover, probMutatie;
int precizie, dimPopulatie, nrEtape, nrBiti;
double fitnessTotal = 0;
vector<double> Interval;
int indice = 0;

struct Individ
{
    string cromozom;
    double x;
    double fitness;
    double ps;

    Individ(string cromozom, double x, double fitness)
        : cromozom(move(cromozom)), x(x), fitness(fitness), ps(0) {}
};

vector<Individ> popInitiala;
vector<Individ> popCurenta;
vector<pair<Individ, int>> popCrossover;

double functie(double x)
{
    return a * x * x + b * x + c;
}

double generareX()
{
    // generam x in intervalul dat
    double x = limInf + static_cast<double>(rand()) / RAND_MAX * (limSup - limInf);
    // il ajustam astfel incat sa aiba nr corect de zecimale
    // ex: 1.2345, precizie 3 => 1234.5 => 1234 => 1.234
    x = round(x * pow(10, precizie)) / pow(10, precizie);
    return x;
}

string codificare(double x)
{
    // d = pasul de discretizare
    double d = (limSup - limInf) / (pow(2, nrBiti) - 1);
    // calculam indexul intervalului de discretizare in care se afla x
    int ind = static_cast<int>((x - limInf) / d);
    // transformam indexul respectiv in binar
    bitset<64> b(ind);
    // substr are ca parametri indicele primului caracter si numarul de caractere din string
    return b.to_string().substr(64 - nrBiti, nrBiti);
}

double decodificare(string cromozom)
{

    bitset<64> bits(cromozom);
    // transformam x in intreg pe baza scrierii in binar
    double x = bits.to_ulong();
    double d = (limSup - limInf) / (pow(2, nrBiti) - 1);
    // afisam capatul din stanga al intervalului cu indexul x
    x = x * d + limInf;
    return x;
}

void generare()
{
    for (int i = 0; i < dimPopulatie; ++i)
    {
        double x = generareX();
        // transformam x in binar
        string cromozom = codificare(x);
        double f = functie(x);
        popInitiala.emplace_back(cromozom, x, f);
        // pastram pentru probabilitatea de selectie
        fitnessTotal += f;
    }

    // pentru fiecare individ, calculam probabilitatea de selectie
    for (auto &individ : popInitiala)
        individ.ps = individ.fitness / fitnessTotal;

    // afisare
    fout << "Populatia initiala\n";
    for (int i = 0; i < popInitiala.size(); ++i)
    {
        fout << i + 1 << ": " << popInitiala[i].cromozom << " x= "<< fixed << setprecision(precizie) << popInitiala[i].x << " f=" << popInitiala[i].fitness << endl;
    }

    fout << "\nProbabilitati selectie\n";
    for (int i = 0; i < popInitiala.size(); ++i)
    {
        fout << "cromozom " << i + 1 << " probabilitate " << popInitiala[i].ps << endl;
    }
}

void intervale()
{
    double sumaPs = 0;
    Interval.push_back(0);
    for (auto &individ : popInitiala)
    {
        //pastram sumele curente pentru a afisa ulterior intervalele
        sumaPs += individ.ps;
        Interval.push_back(sumaPs);
    }
}

void afisareIntervale()
{
    fout << "\nIntervale probabilitati selectie\n";
    for (const auto &suma : Interval)
        fout << suma << " ";
    fout << endl << endl;
}

void popSelectata()
{
    for (int i = 0; i < dimPopulatie; ++i)
    {
        // generam un u de la 0 la 1
        double u = static_cast<double>(rand()) / RAND_MAX;
        // vedem in ce interval apartine u si alegem cromozomul corespunzator pentru selectie
        for (int k = 0; k < Interval.size(); ++k)
        {
            if (u >= Interval[k] && u <= Interval[k+1])
            {
                //individul cu indicele k+1 este selectat (aici avem indicele k pentru ca popInit incepe de la 0)
                popCurenta.push_back(popInitiala[k]);
                if (indice == 0)
                    fout << "u=" << u << fixed << setprecision(precizie) << "  selectam cromozomul " << k + 1 << endl;
                break;
            }
        }
    }
    if (indice == 0)
        fout << endl;
}

void populatieCrossover()
{
    if (indice == 0)
    {
        fout << endl << endl;
        fout << "Probabilitatea de incrucisare " << probCrossover << endl;
        fout << endl;
    }

    for (size_t i = 0; i < popCurenta.size(); ++i)
    {
        // generam un u intre 0 si 1 pentru fiecare cromozom
        double u = static_cast<double>(rand()) / RAND_MAX;

        if (indice == 0)
            fout << i + 1 << ": " << popCurenta[i].cromozom << " u=" << u;

        if (u < probCrossover)
        {
            if (indice == 0)
                fout << "<" << probCrossover << " participa ";
            // am pastrat si indicele fiecarui cromozom pentru a-l afisa corect ulterior
            popCrossover.push_back(make_pair(popCurenta[i], i));
        }
        if (indice == 0)
            fout << endl;
    }
}

void crossover()
{
    if (indice == 0)
        fout << "\nProces de recombinare:\n";
    // se iau perechi de cromozomi
    for (size_t i = 0; i < popCrossover.size(); i += 2)
    {
        //daca avem nr impar de cromozomi, ultimul ramane asa
        if (i + 1 < popCrossover.size())
        {
            // punctul de rupere este un numar random intre 0 si numarul de biti
            int punctRupere = rand() % nrBiti;
            // concatenam prima parte a primului cromozom cu cea de-a doua parte a celuilalt si invers
            string cromInitial1 = popCrossover[i].first.cromozom;
            string cromInitial2 = popCrossover[i + 1].first.cromozom;

            string cromNou1 = cromInitial1.substr(0, punctRupere) + cromInitial2.substr(punctRupere);

            string cromNou2 = cromInitial2.substr(0, punctRupere) + cromInitial1.substr(punctRupere);


            // pastram index pentru afisare
            int index1 = popCrossover[i].second;
            int index2 = popCrossover[i + 1].second;

            // inlocuim cromozomii initiali cu cei noi
            popCurenta[index1].cromozom = cromNou1;
            popCurenta[index2].cromozom = cromNou2;

            // recalculam x si fitness
            popCurenta[index1].x = decodificare(cromNou1);
            popCurenta[index1].fitness = functie(popCurenta[index1].x);

            popCurenta[index2].x = decodificare(cromNou2);
            popCurenta[index2].fitness = functie(popCurenta[index2].x);

            if (indice == 0)
            {
                fout << "Recombinare dintre cromozomul " << index1 + 1 << " cu cromozomul " << index2 + 1 << ":\n";
                fout << cromInitial1 << " " << cromInitial2 << " punct  " << punctRupere << "\n";
                fout << "Rezultat    " << cromNou1 << " " << cromNou2<< "\n\n";
            }
        }
    }
}

void flipBit(string &cromozom, int j) 
{
    if(cromozom[j] == '1') cromozom[j] = '0';
    else cromozom[j] = '1';
}

void afisarePopulatie()
{
    for (int i = 0; i < popCurenta.size(); ++i)
    {
        fout << i + 1 << ": " << popCurenta[i].cromozom
             << " x= " << fixed << setprecision(precizie) << popCurenta[i].x
             << " f=" << popCurenta[i].fitness << endl;
    }
}

void mutatie()
{
    if (indice == 0)
    {
        fout << endl;
        fout << "Probabilitate de mutatie pentru fiecare gena " << probMutatie << "\nAu fost modificati cromozomii:\n";
    }
    //pentru afisarea indicilor
    vector<int> indMod;

    // se ia fiecare individ
    for (int i = 0; i < popCurenta.size(); ++i)
    {
        bool modificat = false;
        // se ia fiecare bit al scrierii in binar
        for (int j = 0; j < nrBiti; ++j)
        {
            // se genereaza o valoare intre 0 si 1; daca e < probabilitate, bit-ul e selectat
            if (static_cast<double>(rand()) / RAND_MAX < probMutatie)
            {
                flipBit(popCurenta[i].cromozom, j);
                modificat = true;
            }
        }

        if (modificat)
        {
            // se retin indicii cromozomilor care au fost modificati
            indMod.push_back(i + 1); 
            // pentru cromozomii modificati, se recalculeaza x si fitness
            popCurenta[i].x = decodificare(popCurenta[i].cromozom);
            popCurenta[i].fitness = functie(popCurenta[i].x);
        }
    }

    if (indice == 0)
        for (int index : indMod)
            fout << index << endl;
}

void elitism()
{
    // inlocuim individul cu cel mai prost fitness din populatia curenta cu individul cu cel mai bun index
    // din populatia trecuta
    Individ individMax = popInitiala[1];
    for (int i = 1; i < dimPopulatie; ++i)
        if (popInitiala[i].fitness > individMax.fitness)
            individMax = popInitiala[i];

    int indFitMin = 0;
    for (int i = 1; i < dimPopulatie; ++i)
        if (popCurenta[i].fitness < popCurenta[indFitMin].fitness)
            indFitMin = i;

    popCurenta[indFitMin] = individMax;
}

double fitnessMaxim()
{
    double fitness_maxim = INT_MIN;
    for (int i = 1; i < dimPopulatie; ++i)
    {
        if (popCurenta[i].fitness > fitness_maxim)
            fitness_maxim = popCurenta[i].fitness;
    }
    return fitness_maxim;
}

int main()
{
    fin >> a >> b >> c;
    fin >> limInf >> limSup;
    fin >> precizie;
    fin >> dimPopulatie;
    fin >> probCrossover;
    fin >> probMutatie;
    fin >> nrEtape;

    // formula pentru a vedea de cati biti avem nevoie pentru a transforma nr scris in zecimal
    nrBiti = ceil(log2((limSup - limInf) * pow(10, precizie)));

    srand((time(0)));
    generare();
    intervale();
    afisareIntervale();
    popSelectata();
    fout << "Dupa selectie:\n";
    afisarePopulatie();
    populatieCrossover();
    crossover();
    fout << "Dupa crossover:\n";
    afisarePopulatie();
    mutatie();
    elitism();
    afisarePopulatie();

    fout << endl
         << "Evolutia maximului:" << endl;
    fout << fitnessMaxim() << endl; // maximul din generatia 1
    for (indice = 1; indice < nrEtape; indice++)
    {
        // reluam procesul, considerand populatia curenta cea initiala pentru urmataorea generatie
        popInitiala = popCurenta;
        popCurenta.clear();
        popCrossover.clear();
        Interval.clear();

        intervale();
        popSelectata();
        populatieCrossover();
        crossover();
        mutatie();
        elitism();
        fout << fitnessMaxim() << endl;
    }
    return 0;
}
