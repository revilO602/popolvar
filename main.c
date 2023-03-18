// Oliver Leontiev
// DSA Zadanie 3, Popolvar
// 2020

//#define CRT_NO_DEPRECATE
//#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct vrchol {
	unsigned int halda_index;
	unsigned char x;
	unsigned char y;
	unsigned int teren : 3;
	struct vrchol* pred; // predosly vrchol, pre vypis od konca
	unsigned int cesta;
}VRCHOL;

typedef struct min_halda {
	VRCHOL** pole;
	int dlzka;
}MIN_HALDA;

typedef struct graf {
	// 2D reprezentacia mapy s ukazovatelmi na strukturi vrcholu na mieste policok. Skala = NULL
	VRCHOL*** vrcholy;
	int vyska;
	int sirka;
	int pocet_princezien;
	VRCHOL* drak;
}GRAF;

MIN_HALDA* vytvorHaldu() {
	MIN_HALDA* halda = (MIN_HALDA*)malloc(sizeof(MIN_HALDA));
	halda->pole = (VRCHOL * *)malloc(sizeof(VRCHOL*));
	halda->dlzka = 0;
	return halda;
}
//vymeni adresy vrcholov v halde ale zachova halda_index
void vymenVrcholi(VRCHOL** a, VRCHOL** b) {
	VRCHOL* tmp = *a;
	unsigned int tmp_i = (*a)->halda_index;
	*a = *b;
	*b = tmp;
	(*b)->halda_index = (*a)->halda_index;
	(*a)->halda_index = tmp_i;
}
// opravi haldu smerom hore od indexu, pouzite pri vlozeni vrcholu na koniec haldy
void heapify_up(MIN_HALDA* halda, int index) {
	int i = index;
	while (i > 0 && (halda->pole[(i - 1) / 2])->cesta > (halda->pole[i])->cesta) {
		vymenVrcholi(&(halda->pole)[i], &(halda->pole[(i - 1) / 2]));
		i = (i - 1) / 2;
	}
}
// vlozi na koniec haldy a prenesie smerom hore na spravne miesto
void vlozHalda(MIN_HALDA* halda, VRCHOL* vrchol) {
	halda->dlzka++;
	int i = halda->dlzka - 1;
	halda->pole = realloc(halda->pole, (halda->dlzka) * sizeof(VRCHOL*));
	(halda->pole)[i] = vrchol;
	vrchol->halda_index = i;
	heapify_up(halda, i);
}
// opravi haldu smerom dole, pouzite pri extrakcii z haldy
void heapify_down(MIN_HALDA* halda, int index) {
	int i = index;
	while (2 * i + 1 < halda->dlzka &&
		(halda->pole[i]->cesta > halda->pole[2 * i + 1]->cesta ||
		(2 * i + 2 < halda->dlzka && halda->pole[i]->cesta > halda->pole[2 * i + 2]->cesta))) {
		if (halda->pole[i]->cesta > halda->pole[2 * i + 1]->cesta &&
			(2 * i + 2 >= halda->dlzka || halda->pole[2 * i + 1]->cesta <= halda->pole[2 * i + 2]->cesta)) {
			vymenVrcholi(&(halda->pole[i]), &(halda->pole[2 * i + 1]));
			i = 2 * i + 1;
		}
		else {
			vymenVrcholi(&(halda->pole[i]), &(halda->pole[2 * i + 2]));
			i = 2 * i + 2;
		}
	}
}
// vytiahne zaciatok haldy (minimum), a vymeni ho za koniec. Opravi haldu smerom dole
VRCHOL* vytiahniHalda(MIN_HALDA* halda) {
	VRCHOL* min;
	int i = halda->dlzka - 1;
	if (halda->dlzka == 0)return NULL;
	min = halda->pole[0];
	halda->pole[0] = halda->pole[i];
	halda->pole[0]->halda_index = 0;
	halda->dlzka--;
	halda->pole = realloc(halda->pole, (halda->dlzka) * sizeof(VRCHOL*));
	heapify_down(halda, 0);
	return min;
}

void uvolniHaldu(MIN_HALDA* halda) {
	//strukturi v halde sa uvolnia pri uvolneni grafu
	free(halda->pole);
	free(halda);
}
//zakoduje teren do cisla 0-7 podla policka a poctu princezien
void nastavTeren(VRCHOL* vrchol, unsigned char teren, int pocet_princezien) {
	// drak je 2, princezne su 3-7
	switch (teren) 
	{
	case 'C': vrchol->teren = 0;
		break;
	case 'H': vrchol->teren = 1;
		break;
	case 'D': vrchol->teren = 2;
		break;
	case 'P': vrchol->teren = pocet_princezien+3;
		break;
	default:printf("neznamy symbol na mape\n");
	}
}
//vrati cenu cesty podla terenu
unsigned int zistiCestu(unsigned int teren) {
	// hustina je 1 a iba ta stoji viac
	if (teren == 1)
		return 2;
	return 1;
}

//vytvori a nastavi novy vrchol
VRCHOL* initVrchol(unsigned char x, unsigned char y, char teren,unsigned char start_x, unsigned char start_y, int pocet_princezien) {
	VRCHOL* vrchol = (VRCHOL*)malloc(sizeof(VRCHOL));
	nastavTeren(vrchol,teren,pocet_princezien);
	vrchol->x = x;
	vrchol->y = y;
	
	if (x == start_x && y == start_y) {
		vrchol->cesta = zistiCestu(teren);
		vrchol->pred = vrchol;//ukazuje sam na seba, len aby neukazoval na NULL. NULL znamena ze sa sem neda dostat z ineho vrchola
	}
	else {
		vrchol->cesta = INFINITY;
		vrchol->pred = NULL;
	}
	return vrchol;
}
//vytvori graf z mapy, a kazdy validny vrchol (nie skala) vlozi do haldy
void initGraf(char** mapa,GRAF** graf,unsigned char vyska, unsigned char sirka, MIN_HALDA* halda) {
	unsigned int y, x,p=0,d=0;

	(*graf) = (GRAF*)malloc(sizeof(GRAF));
	(*graf)->vyska = vyska;
	(*graf)->sirka = sirka;
	(*graf)->vrcholy = (VRCHOL * **)malloc(vyska * sizeof(VRCHOL * *));
	for (y = 0; y < vyska; y++) {
		(*graf)->vrcholy[y] = (VRCHOL * *)malloc(sirka * sizeof(VRCHOL*));
		for (x = 0; x < sirka; x++) {
			if (mapa[y][x] != 'N') {
				(*graf)->vrcholy[y][x] = initVrchol(x, y, mapa[y][x], 0, 0,p);
				vlozHalda(halda, (*graf)->vrcholy[y][x]);
				if (mapa[y][x] == 'D') {
					d++;
					(*graf)->drak = (*graf)->vrcholy[y][x];
				}
				if (mapa[y][x] == 'P') p++;
			}
			else (*graf)->vrcholy[y][x] = NULL;
		}
	}
	//v pocet_princezien uchovava aj informaciu a platnosti mapy
	if (d > 1)
		(*graf)->pocet_princezien = -2;
	else if (d < 1)
		(*graf)->pocet_princezien = -1;
	else (*graf)->pocet_princezien = p;
		
}

void uvolniGraf(GRAF* graf) {
	int i, j;

	for (i = 0; i < graf->vyska; i++) {
		for (j = 0; j < graf->sirka; j++) {
			free(graf->vrcholy[i][j]);
		}
		free(graf->vrcholy[i]);
	}
	free(graf->vrcholy);
	free(graf);

}

VRCHOL* dijkstra(MIN_HALDA* halda, GRAF* graf, unsigned int end_terrain);

//obnovi haldu s novym zaciatkom pre dalsieho dijkstru
void resetHalda(MIN_HALDA* halda, GRAF* graf, unsigned char start_x, unsigned char start_y) {
	int i, j;
	free(halda->pole);
	halda->pole = (VRCHOL * *)malloc(sizeof(VRCHOL*));
	halda->dlzka = 0;
	for (i = 0; i < graf->vyska; i++) {
		for (j = 0; j < graf->sirka; j++) {
			if (graf->vrcholy[i][j] != NULL) {
				if (graf->vrcholy[i][j]->x == start_x && graf->vrcholy[i][j]->y == start_y) {
					graf->vrcholy[i][j]->cesta = 0;
					graf->vrcholy[i][j]->pred = graf->vrcholy[i][j];
				}
				else {
					graf->vrcholy[i][j]->cesta = INFINITY;
					graf->vrcholy[i][j]->pred = NULL;
				}
				vlozHalda(halda, graf->vrcholy[i][j]);
			}
		}
	}
}
//spocita cestu od konca po zaciatok po dijkstrovy
void spocitajCestu(VRCHOL* koniec, int* dlzka_cesty, unsigned char start_x, unsigned char start_y) {
	while (!(koniec->x == start_x && koniec->y == start_y)) {
		(*dlzka_cesty)++;
		koniec = koniec->pred;
	}
}

// vyskusa jednu permutaciu princezien a ak je zatial najlacnejsia prepise doterajsiu najlacnejsiu
// princezne=jedna permutacia v poli, cesta=suradnice najlacnejsej cesty v poli, najlacnejsia=cena najlacnejsej cesty
void skusPermutaciu(MIN_HALDA* halda, GRAF* graf, int* princezne, int** cesta, int* najlacnejsia, int* p_dlzka) {
	int i,j=0,k=0,dlzka=0,cena=0;
	VRCHOL* koniec;
	VRCHOL* start = graf->drak;
	int* akt_cesta=NULL;
	for (i = 0; i < graf->pocet_princezien; i++) {
		resetHalda(halda, graf, start->x, start->y);
		koniec = dijkstra(halda, graf, princezne[i]);
		if (koniec == NULL) {
			*p_dlzka = -2; // -2 znamena ze cesta je zablokovana, vsetko sa prerusi 
			free(akt_cesta);
			return;
		}
		cena += koniec->cesta;
		k = dlzka;
		spocitajCestu(koniec, &dlzka, start->x, start->y);
		j = dlzka;
		akt_cesta = realloc(akt_cesta, sizeof(int) * dlzka * 2);
		start = koniec;
		while (j > k) {
			akt_cesta[(--j) * 2] = koniec->x;
			akt_cesta[j * 2 + 1] = koniec->y;
			koniec = koniec->pred;
		}
	}
	// bud uvolni akt_cestu alebo uvolni cestu a nasmeruje smernik na akt_cestu
	if (cena < *najlacnejsia || *najlacnejsia == -1) {// najlacnejsia je -1 na zaciatku
		*najlacnejsia = cena;
		free(*cesta);
		*cesta = akt_cesta;
		*p_dlzka = dlzka;
	}
	else free(akt_cesta);
}
// vymeni integery
void vymen(int* a, int* b) {
	int tmp = *a;
	*a = *b;
	*b = tmp;
}
//vyuzitie Heapovho algoritmu na zistenie vsetkych permutacii, permutacie su z cisel terenu princezien (3-7 pre 5 princezien)
void permutuj(MIN_HALDA* halda, GRAF* graf, int** princezne, int pocet_princezien, int** cesta, int* najlacnejsia, int* p_dlzka) {
	int i;

	if (*p_dlzka == -2) //ak je jedna cesta zablokovana, tak vsetky
		return;
	if (pocet_princezien == 1) {
		skusPermutaciu(halda, graf, *princezne, cesta, najlacnejsia,p_dlzka);
	}
	else {
		permutuj(halda, graf, princezne, pocet_princezien - 1, cesta, najlacnejsia, p_dlzka);
		for (i = 0; i < pocet_princezien - 1; i++) {
			if (pocet_princezien % 2 == 0) {
				vymen(&(*princezne)[i], &(*princezne)[pocet_princezien - 1]);
			}
			else {
				vymen(&(*princezne)[0], &(*princezne)[pocet_princezien - 1]);
			}
			permutuj(halda, graf, princezne, pocet_princezien - 1, cesta, najlacnejsia, p_dlzka);
		}
	}
}

// vytvor pole pre permutacie a vsetky vyskusaj
int vsetkyPermutacie(MIN_HALDA* halda, GRAF* graf, int** cesta) {
	int* princezne = (int*)malloc(sizeof(int)*graf->pocet_princezien);
	int i,p_dlzka=0, najlacnejsia=-1;
	for (i = 0; i < graf->pocet_princezien;i++) {
		princezne[i] = i + 3;
	}
	permutuj(halda, graf, &princezne, graf->pocet_princezien, cesta, &najlacnejsia,&p_dlzka);
	free(princezne);
	return p_dlzka;
}

// dijkstrov algoritmus s vyuzitim min. haldy. koniec_teren je teren policka ku ktoremu hladame cestu
// 2 = drak, 3-7= princezne
VRCHOL* dijkstra(MIN_HALDA* halda,GRAF* graf, unsigned int koniec_teren) {
	VRCHOL* akt;
	VRCHOL* akt_sused=NULL;
	int i;
	
	while (1) {
		akt = vytiahniHalda(halda);
		if (akt == NULL || akt->pred==NULL) {
			printf("Ciel je zablokovany\n");
			return NULL;
		}
		for (i = 0; i < 4; i++) { //vyhodnot cestu do vsetkych susedov
			switch (i) {
			case 0: 
				if (akt->x + 1 < graf->sirka)
					akt_sused = graf->vrcholy[akt->y][akt->x+1];
				else akt_sused = NULL;
				break;
			case 1:
				if (akt->x - 1 >= 0)
					akt_sused = graf->vrcholy[akt->y][akt->x-1];
				else akt_sused = NULL;
				break;
			case 2:
				if (akt->y + 1 < graf->vyska)
					akt_sused = graf->vrcholy[akt->y+1][akt->x];
				else akt_sused = NULL;
				break;
			case 3:
				if (akt->y-1 >= 0) 
					akt_sused = graf->vrcholy[akt->y-1][akt->x];
				else akt_sused = NULL;
				break;
			}
			if (akt_sused == NULL)continue;
			if (akt_sused->cesta > akt->cesta + zistiCestu(akt_sused->teren)) {
				akt_sused->cesta = zistiCestu(akt_sused->teren) + akt->cesta;
				akt_sused->pred = akt;
			}
			heapify_up(halda, akt_sused->halda_index);
		}
		if (akt->teren == koniec_teren) break;
	}
	return akt;
}

int skontrolujMapu(int p) {
	if (p > 5) {
		printf("Na mape je privela princezien\n");
		return 0;
	}
	if (p ==-1) {
		printf("Na mape nie je drak\n");
		return 0;
	}
	if (p == -2) {
		printf("Na mape je viac ako jeden drak\n");
		return 0;
	}
	if (p == 0) {
		printf("Na mape nie su princezne\n");
		return 0;
	}
	return 1;
}
int* zabiDraka(GRAF* graf, MIN_HALDA* halda,int* dlzka_cesty) {
	VRCHOL* koniec;
	int* cesta;
	int i;
	if ((koniec = dijkstra(halda, graf, 2)) == NULL) {
		*dlzka_cesty = 0;
		return NULL;
	}
	(*dlzka_cesty)++;
	spocitajCestu(koniec, dlzka_cesty, 0, 0);
	cesta = (int*)malloc(sizeof(int) * *dlzka_cesty * 2);
	i = *dlzka_cesty;
	cesta[(--i) * 2] = koniec->x;
	cesta[i * 2 + 1] = koniec->y;
	while (i > 0) {
		koniec = koniec->pred;
		cesta[(--i) * 2] = koniec->x;
		cesta[i * 2 + 1] = koniec->y;
	}
	return cesta;
}

int* zachran_princezne(char** mapa, int n, int m, int t, int* dlzka_cesty) {
	GRAF* graf;
	int* cesta=NULL,*p_cesta=NULL;
	int i,k,p_dlzka=0;
	MIN_HALDA* halda = vytvorHaldu();

	initGraf(mapa, &graf, n, m, halda);//vytvori graf a nacita vrcholy do haldy
	if (skontrolujMapu(graf->pocet_princezien) == 0)return NULL;
	cesta = zabiDraka(graf, halda, dlzka_cesty);
	if (cesta == NULL) {
		printf("Neda sa dostat ku drakovi\n");
		return NULL;
	}
	p_dlzka = vsetkyPermutacie(halda,graf,&p_cesta);
	if (p_dlzka == -2) {
		free(cesta);
		*dlzka_cesty = 0;
		printf("Neda sa dostat ku princeznej\n");
		return NULL;
	}
	k = *dlzka_cesty*2;
	*dlzka_cesty += p_dlzka;
	cesta = realloc(cesta, sizeof(int) * *dlzka_cesty * 2);
	for (i = 0; i < p_dlzka *2; i++) {
		cesta[k + i] = p_cesta[i];
	}
	free(p_cesta);
	uvolniGraf(graf);
	uvolniHaldu(halda);
	return cesta;
}

int main()
{
	char** mapa;
	int i, test, dlzka_cesty, cas, * cesta;
	int n = 0, m = 0;
	unsigned int t = 0;
	FILE* f, *f_moj=NULL;
	while (1) {
		printf("Zadajte cislo testu (0 ukonci program):\n");
		scanf("%d", &test);
		dlzka_cesty = 0;
		n = m = t = 0;
		switch (test) {
		case 0://ukonci program
			if (f_moj)
				fclose(f_moj);
			return 0;
		case 1://nacitanie mapy zo suboru
			f = fopen("test.txt", "r");
			if (f)
				fscanf(f, "%d %d %d", &n, &m, &t);
			else
				continue;
			mapa = (char**)malloc(n * sizeof(char*));
			for (i = 0; i < n; i++) {
				mapa[i] = (char*)malloc(m * sizeof(char));
				for (int j = 0; j < m; j++) {
					char policko = fgetc(f);
					if (policko == '\n') policko = fgetc(f);
					mapa[i][j] = policko;
				}
			}
			fclose(f);
			cesta = zachran_princezne(mapa, n, m, t, &dlzka_cesty);
			break;
		case 2://nacitanie preddefinovanej mapy
			n = 10;
			m = 10;
			t = 12;
			mapa = (char**)malloc(n * sizeof(char*));
			mapa[0] = "CCHCNHCCHN";
			mapa[1] = "NNCCCHHCCC";
			mapa[2] = "DNCCNNHHHC";
			mapa[3] = "CHHHCCCCCC";
			mapa[4] = "CCCCCNHHHH";
			mapa[5] = "PCHCCCNNNN";
			mapa[6] = "NNNNNHCCCC";
			mapa[7] = "CCCCCPCCCC";
			mapa[8] = "CCCNNHHHHH";
			mapa[9] = "HHHPCCCCCC";
			cesta = zachran_princezne(mapa, n, m, t, &dlzka_cesty);
			break;
		case 3: 
			if (!f_moj) {
				f_moj = fopen("moj_test.txt", "r");
			}
			if (f_moj)
				fscanf(f_moj, "%d %d %d", &n, &m, &t);
			else
				continue;
			mapa = (char**)malloc(n * sizeof(char*));
			for (i = 0; i < n; i++) {
				mapa[i] = (char*)malloc(m * sizeof(char));
				for (int j = 0; j < m; j++) {
					char policko = fgetc(f_moj);
					if (policko == '\n') policko = fgetc(f_moj);
					mapa[i][j] = policko;
				}
			}
			cesta = zachran_princezne(mapa, n, m, t, &dlzka_cesty);
			break;
		default:
			continue;
		}
		cas = 0;
		for (i = 0; i < dlzka_cesty; i++) {
			printf("%d %d\n", cesta[i * 2], cesta[i * 2 + 1]);
			if (mapa[cesta[i * 2 + 1]][cesta[i * 2]] == 'H')
				cas += 2;
			else
				cas += 1;
			if (mapa[cesta[i * 2 + 1]][cesta[i * 2]] == 'D') {
				if (cas > t) {
					printf("Nestihol si zabit draka!\n");
					free(cesta);
					if (test != 2)
						for (i = 0; i < n; i++) {
							free(mapa[i]);
						}
					free(mapa);
					return 0;
				}
				else {
					t = INFINITY; // ak zabil draka uz ma nekonecno casu
				}
			}
				
			if (mapa[cesta[i * 2 + 1]][cesta[i * 2]] == 'N')
				printf("Prechod cez nepriechodnu prekazku!\n");
			if (i > 0 && abs(cesta[i * 2 + 1] - cesta[(i - 1) * 2 + 1]) + abs(cesta[i * 2] - cesta[(i - 1) * 2]) > 1)
				printf("Neplatny posun Popolvara!\n");
		}
		printf("%d\n", cas);
		free(cesta);
		if (test != 2)
			for (i = 0; i < n; i++) {
				free(mapa[i]);
			}
		free(mapa);
	}
	if (f_moj)
		fclose(f_moj);
	return 0;
}