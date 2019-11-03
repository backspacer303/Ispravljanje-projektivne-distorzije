#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

#define cimg_use_jpeg
#define cimg_use_png
#include "eigen/Eigen/Eigen"
#include "eigen/Eigen/SVD"
#include "CImg-2.7.4/CImg.h"


std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> 
unesi_tacke(void){
    
    double x;
    std::vector<double> tmp;
    
    std::vector<std::vector<double>> tacke;
    std::vector<std::vector<double>> slike;
    
    int br_tacaka;
    
    std::cout << "Uneti broj tacaka: " << std::endl;
    std::cin >> br_tacaka;
    
    int i;
    std::cout << "Unos tacaka:" << std::endl;
    for(i = 1; i<=br_tacaka*3; i++){
        
        std::cin >> x;
        tmp.push_back(x);
        
        if(i%3 == 0){
            tacke.push_back(tmp);
            tmp.clear();
        }
    }
    
    std::cout << "Unos slika:" << std::endl;
    for(i = 1; i<=br_tacaka*3; i++){
        
        std::cin >> x;
        tmp.push_back(x);
        
        if(i%3 == 0){
            slike.push_back(tmp);
            tmp.clear();
        }
    }
    
    return std::make_pair(tacke, slike);
}

Eigen::MatrixXd naivni(void){
    
    double x;
    std::vector<double> tmp;
    
    std::vector<std::vector<double>> tacke;
    
    std::cout << "Unos tacaka:" << std::endl;
    for(int i=1; i<=12; i++){
        
        std::cin >> x;
        tmp.push_back(x);
        
        if(i%3 == 0){
            tacke.push_back(tmp);
            tmp.clear();            
        }        
    }
    
   
    Eigen::MatrixXd delta;
    delta.resize(3,3);
    delta << tacke[0][0], tacke[1][0], tacke[2][0],
             tacke[0][1], tacke[1][1], tacke[2][1],
             tacke[0][2], tacke[1][2], tacke[2][2];
          
    Eigen::MatrixXd delta1;
    delta1.resize(3,3);
    delta1 << tacke[3][0], tacke[1][0], tacke[2][0],
              tacke[3][1], tacke[1][1], tacke[2][1],
              tacke[3][2], tacke[1][2], tacke[2][2];
  
    Eigen::MatrixXd delta2;
    delta2.resize(3,3);    
    delta2 << tacke[0][0], tacke[3][0], tacke[2][0],
              tacke[0][1], tacke[3][1], tacke[2][1],
              tacke[0][2], tacke[3][2], tacke[2][2];
   
    Eigen::MatrixXd delta3;
    delta3.resize(3,3);    
    delta3 << tacke[0][0], tacke[1][0], tacke[3][0],
              tacke[0][1], tacke[1][1], tacke[3][1],
              tacke[0][2], tacke[1][2], tacke[3][2];
             
    //determinante         
    double d  = delta.determinant();         
    double d1 = delta1.determinant();
    double d2 = delta2.determinant();
    double d3 = delta3.determinant();
    
    //lambde
    double l1 = d1/d;
    double l2 = d2/d;
    double l3 = d3/d;
    
    
    Eigen::MatrixXd p1;
    p1.resize(3,3);
        
    p1 << tacke[0][0]*l1, tacke[1][0]*l2, tacke[2][0]*l3, 
          tacke[0][1]*l1, tacke[1][1]*l2, tacke[2][1]*l3,
          tacke[0][2]*l1, tacke[1][2]*l2, tacke[2][2]*l3;

    Eigen::MatrixXd p1_i = p1.inverse();
  
    
    std::vector<std::vector<double>> tacke2;
    
    std::cout << "Unos slika tacaka: " << std::endl;    
    for(int i=1; i<=12; i++){
    
        std::cin >> x;
        tmp.push_back(x);
        
        if(i%3 == 0){
            tacke2.push_back(tmp);
            tmp.clear();            
        }        
    }
    
    delta << tacke2[0][0], tacke2[1][0], tacke2[2][0],
             tacke2[0][1], tacke2[1][1], tacke2[2][1],
             tacke2[0][2], tacke2[1][2], tacke2[2][2];
             
    delta1 << tacke2[3][0], tacke2[1][0], tacke2[2][0],
              tacke2[3][1], tacke2[1][1], tacke2[2][1],
              tacke2[3][2], tacke2[1][2], tacke2[2][2];
              
    delta2 << tacke2[0][0], tacke2[3][0], tacke2[2][0],
              tacke2[0][1], tacke2[3][1], tacke2[2][1],
              tacke2[0][2], tacke2[3][2], tacke2[2][2];
              
    delta3 << tacke2[0][0], tacke2[1][0], tacke2[3][0],
              tacke2[0][1], tacke2[1][1], tacke2[3][1],
              tacke2[0][2], tacke2[1][2], tacke2[3][2];
    
    d  = delta.determinant();
    d1 = delta1.determinant();
    d2 = delta2.determinant();
    d3 = delta3.determinant();
    
    l1 = d1/d;
    l2 = d2/d;
    l3 = d3/d;
    
    
    
    Eigen::MatrixXd p2;
    p2.resize(3,3);
        
    p2 << tacke2[0][0]*l1, tacke2[1][0]*l2, tacke2[2][0]*l3, 
          tacke2[0][1]*l1, tacke2[1][1]*l2, tacke2[2][1]*l3,
          tacke2[0][2]*l1, tacke2[1][2]*l2, tacke2[2][2]*l3;
    
    Eigen::MatrixXd g;
    g = p2*p1_i;
    
    return g;
}

Eigen::MatrixXd DLT(std::vector<std::vector<double>> tacke, std::vector<std::vector<double>> slike){

    int br_tacaka = tacke.size();
    
    Eigen::MatrixXd A;
    A.resize(2*br_tacaka, 9);
    

    for(int i=0; i<br_tacaka; i++){
        
        A.row(i*2) << 0, 0, 0, 
                     (-1)*slike[i][2]*tacke[i][0], (-1)*slike[i][2]*tacke[i][1], (-1)*slike[i][2]*tacke[i][2],
                     slike[i][1]*tacke[i][0], slike[i][1]*tacke[i][1], slike[i][1]*tacke[i][2];            
            
        
        A.row(i*2+1) << slike[i][2]*tacke[i][0], slike[i][2]*tacke[i][1], slike[i][2]*tacke[i][2],
                        0, 0, 0,
                        (-1)*slike[i][0]*tacke[i][0], (-1)*slike[i][0]*tacke[i][1], (-1)*slike[i][0]*tacke[i][2];
    }
    
    
    std::cout << std::endl << "Matrica A: \n" << A << std::endl << std::endl;
    
    
    
    Eigen::JacobiSVD<Eigen::MatrixXd> s(A, Eigen::ComputeFullV | Eigen::ComputeFullU );     
    //s.computeV();
    Eigen::MatrixXd V = s.matrixV();
    std::cout << "Matrica V: \n" << V << std::endl << std::endl;
    
    
    Eigen::MatrixXd P;
    P.resize(3,3);
    P << V(0,8), V(1,8), V(2,8),
         V(3,8), V(4,8), V(5,8),
         V(6,8), V(7,8), V(8,8);
    
    
    
    return P;
}

std::vector<std::vector<double>> toAfine(std::vector<std::vector<double>> t){
    
    int n = t.size();
    
    std::vector<std::vector<double>> tacke_a;
    std::vector<double> tmp;
    
    for(int i=0; i<n; i++){
        
        tmp.push_back(t[i][0]/t[i][2]);
        tmp.push_back(t[i][1]/t[i][2]);
        
        tacke_a.push_back(tmp);
        tmp.clear();
    }
    return tacke_a;
}
std::vector<std::vector<double>> toHomogeneous(std::vector<std::vector<double>> t){
    
    int n = t.size();
    
    std::vector<std::vector<double>> tacke_h;
    std::vector<double> tmp;
    
    for(int i=0; i<n; i++){
        
        tmp.push_back(t[i][0]);
        tmp.push_back(t[i][1]);
        tmp.push_back(1);
        
        tacke_h.push_back(tmp);
        tmp.clear();
    }
    
    return tacke_h;
    
}

std::pair<Eigen::MatrixXd, std::vector<std::vector<double>>> 
    normalizuj(std::vector<std::vector<double>> tacke){
    
        
    tacke = toAfine(tacke);
    
    int n = tacke.size();
    
    std::vector<double> C;
    C.resize(2);
    
    for(int i = 0; i<n; i++){
        C[0] += tacke[i][0]; 
        C[1] += tacke[i][1];
        
    }
    
    C[0] /= n;
    C[1] /= n;
    
    
    Eigen::MatrixXd G;
    G.resize(3,3);
    
    G << 1, 0, -C[0],
         0, 1, -C[1],
         0, 0, 1;
         
    std::vector<std::vector<double>> tacke_t;
    std::vector<double> tmp;
    double x;
    
    //translacija tacaka
    for(int i=0; i<n; i++){
       
        tmp.push_back(tacke[i][0]+G(0,2));
        tmp.push_back(tacke[i][1]+G(1,2));
                
        tacke_t.push_back(tmp);
        tmp.clear();
    }
    
    //racunanje lambde
    double lamda = 0;
    for(int i=0; i<n; i++){
        lamda += sqrt( pow(tacke_t[i][0],2) + pow(tacke_t[i][1],2) ); 
    }
    
    lamda /= n;
    
    
    //racunanje matrice S
        
    double k = sqrt(2)/lamda;
    
    Eigen::MatrixXd S;
    S.resize(3,3);
    
    S << k, 0, 0,
         0, k, 0,
         0, 0, 1;
    
    //racunanje matrice T
    Eigen::MatrixXd T;
    T.resize(3,3);
    T = S*G;
         
    //homotetija tacaka
    std::vector<std::vector<double>> tacke_t_h;
    tmp.clear();
    for(int i=0; i<n; i++){
        
        tmp.push_back(tacke_t[i][0]*k);
        tmp.push_back(tacke_t[i][1]*k);
        
        tacke_t_h.push_back(tmp);
        tmp.clear();
    }
    
    std::vector<std::vector<double>> tacke_th_H;
    tacke_th_H = toHomogeneous(tacke_t_h);
    
    
    return std::make_pair(T, tacke_th_H);     
}

Eigen::MatrixXd 
DLT_norm(std::vector<std::vector<double>> tacke, std::vector<std::vector<double>> slike){
    
    
    
    
    auto N1 = normalizuj(tacke);  //NORMALIZACIJA ORIGINALA
    auto N2 = normalizuj(slike);  //NORMALIZACIJA SLIKA
    
    Eigen::MatrixXd P1;     
    Eigen::MatrixXd T1;
    Eigen::MatrixXd T2;
    Eigen::MatrixXd P;
    
    P1 = DLT(N1.second, N2.second);
    
    T1 = N1.first;
    T2 = N2.first.inverse();
    
    P = T2*P1*T1;
    
    
    std::cout << std::endl << "Matrica T1:\n" << T1 << std::endl << std::endl;
    std::cout << std::endl << "Matrica T2:\n" << N2.first << std::endl << std::endl;
    std::cout << std::endl << "Matrica P1:\n" << P1 << std::endl << std::endl;
    std::cout << std::endl << "Matrica P = T2_inv*P1*T1:\n" << P << std::endl << std::endl;
    
    return P;
}


void ispravi_distorziju(){
    
    int odabir_sike;
    
    while(1){
        
        
        std::cout << std::endl << "Odaberite sliku (1-5): ";
        std::cin >> odabir_sike;
        
        if(odabir_sike == 1 || odabir_sike == 2 || odabir_sike == 3 || odabir_sike == 4 || odabir_sike == 5 )
            break;
        else
            std::cout << std::endl << std::endl << "Greska! Pokusajte ponovo." << std::endl;
    }
    
   
   cimg_library::CImg<unsigned char> img;
    
    switch(odabir_sike){
        
        case 1:
            img.assign("ulazne_slike/zgrada2.jpg");
            break;
            
        case 2:
            img.assign("ulazne_slike/zgrada1.png");
            break;
            
        case 3:
            img.assign("ulazne_slike/Two-Points-Perspective-3.jpg");
            break;
            
        case 4:
            img.assign("ulazne_slike/crtez1.jpg");
            break;
            
        case 5:
            img.assign("ulazne_slike/crtez2.jpg");
            break; 
        
    }
    
        
    
        
    std::cout << img.width() << "..x.." <<img.height() << std::endl;
    
    float w = img.width();
    float h = img.height();
    
    double ratio = floor(w/h);
    
    //cetvrtu sliku necemo da resize-ujemo: nece lepo da je skalira iz nekog razloga
    if(odabir_sike != 4)
        img.resize(900, 900/ratio);
    
    
    cimg_library::CImgDisplay display1 (img, "Odaberite tacke klikom na sliku i pritisnite ENTER!");
    
    
    
    std::vector<std::vector<double>> odabrane_tacke;
    std::vector<double> tmp;
    
    
    
    
    std::cout << "Koordinate odabranih tacaka:" <<std::endl;
    while(!display1.is_closed() &&  !display1.is_keyENTER()){
        
        
        display1.wait();
        
        
        if(display1.button() && display1.mouse_y() >= 0 && display1.mouse_x() >= 0){
            
            const int x = display1.mouse_x();
            const int y = display1.mouse_y();
            
            tmp.push_back(x);
            tmp.push_back(y);
            odabrane_tacke.push_back(tmp);
            tmp.clear();
            
            std::cout << x << " " << y << std::endl;
            
        }
        
       
        
    }

    display1.close();
    
    
    /*
    std::cout << "++++++++++++++++++++++\n";
    
    for(int i = 0; i<4; i++){
        
        for(int j =0; j<2; j++){
            std::cout << odabrane_tacke[i][j]  << " ";
        }
        
        std::cout << std::endl;
    }
    
    std::cout << "++++++++++++++++++++++\n";
    */
    
    std::vector<std::vector<double>> odabrane_tacke_H;
    odabrane_tacke_H = toHomogeneous(odabrane_tacke);
    
    
    
    std::vector<std::vector<double>> slike_tacaka;
    int br_tacaka = odabrane_tacke.size();
    
    slike_tacaka.resize(br_tacaka);
    for(int i=0; i<br_tacaka; i++)
        slike_tacaka[i].resize(3);    
    
    
    std::cout << "Uneti koordinate slika:" << std::endl;
    for(int i=0; i<br_tacaka; i++){
        
        std::cin >> slike_tacaka[i][0];
        std::cin >> slike_tacaka[i][1];        
    }
    
    /*
    std::cout << ">>>>>>>>>>>>>>>>>>>\n";
    
    for(int i = 0; i<4; i++){
        
        for(int j =0; j<2; j++){
            std::cout << slike_tacaka[i][j]  << " ";
        }
        
        std::cout << std::endl;
    }        
    std::cout << ">>>>>>>>>>>>>>>>>>>\n";
    */
    

    
    std::vector<std::vector<double>> slike_tacaka_H;
    slike_tacaka_H = toHomogeneous(slike_tacaka);
    
    
    
    
    Eigen::MatrixXd P;
    P.resize(3,3);
    
    P = DLT_norm(odabrane_tacke_H, slike_tacaka_H);
    
    Eigen::MatrixXd Pinv;
    Pinv.resize(3,3);
    Pinv = P.inverse();
    
    
    
    int width  = img.width();
    int height = img.height();    
    cimg_library::CImg<unsigned char> img2 (width, height, 1, 3, 0);    
    
    
    
    Eigen::MatrixXd pixel;
    pixel.resize(3,1);

    
    Eigen::MatrixXd px_slika;
    px_slika.resize(3,1);
    
    int xo, yo;
    
    for(int x=0; x<width; x++){
        for(int y=0; y<height; y++){
            
            pixel << x,y,1;
            
            px_slika = Pinv*pixel;
            
            
            if(px_slika(2,0) != 0){
                xo = floor(px_slika(0,0)/px_slika(2,0));
                yo = floor(px_slika(1,0)/px_slika(2,0));
            }
            else
                continue;
            
            
            
            if( xo < 0 || xo > width || yo < 0 || yo > height )
                continue;
            else{
                
                img2(x,y,0,0) = (int)img(xo, yo, 0, 0);
                img2(x,y,0,1) = (int)img(xo, yo, 0, 1);
                img2(x,y,0,2) = (int)img(xo, yo, 0, 2);
            }
            
        }
    }
    
    
    
    
    
    
    
    
    
    cimg_library::CImgList<unsigned char> lista_slika (img, img2);
    
    cimg_library::CImgDisplay disp2 (lista_slika, "Otklonjena distorzija :D");
    
    disp2.move(-50, -50);
    
    while(!disp2.is_closed() &&  !disp2.is_keyESC()){
                
        display1.wait(1);        
                        
    }
    
       
}



int main(){
      
    while(1){
        
        int odgovor;
        std::cout << "Izaberite opciju: \n1 - testiranje algoritama\n2 - ispravljanje distorzije \nOdgovor: ";
        std::cin >> odgovor;
        std::cout << std::endl << std::endl;
        
        if(odgovor == 1){
        
            std::cout << "******************* NAIVNI **********************" << std::endl << std::endl;
            
            Eigen::MatrixXd P = naivni();
            std::cout << "Matrica projekcije (naivni): \n" << P << std::endl;
            
            std::cout << "*************************************************" << std::endl << std::endl;
            
                    
            std::cout << "******************** DLT ************************" << std::endl << std::endl;
            
            auto u = unesi_tacke();
            Eigen::MatrixXd D = DLT(u.first, u.second);
            
            std::cout << "Matrica projekcije (naivni): \n" << P << std::endl << std::endl;
            std::cout << "Matrica projekcije (DLT): \n" << D << std::endl << std::endl;
            
            //WARNING: Ako je D(0,0) = 0 javice gresku zbog deljenja nulom
            Eigen::MatrixXd PR = (D/D(0,0))*P(0,0);
            std::cout << "Provera (DLT->naivni): \n" << PR << std::endl;
            
            std::cout << std::endl <<"*************************************************" << std::endl << std::endl;
            

            
            std::cout << "******************** DLT_NORM ************************" << std::endl << std::endl;
            
            auto unos = unesi_tacke();
            Eigen::MatrixXd M = DLT_norm(unos.first, unos.second);    
            Eigen::MatrixXd PM = (M/M(0,0))*P(0,0);    
            std::cout << "Matrica projekcije (naivni): \n" << P << std::endl << std::endl;    
            std::cout << "Matrica projekcije (DLT_norm): \n" << M << std::endl << std::endl;
            std::cout << "Provera (DLT_norm->naivni): \n" << PM << std::endl << std::endl;
            
            std::cout << "******************************************************" << std::endl << std::endl;
            
            break;
        }
        
        else if(odgovor == 2){
            
            ispravi_distorziju();
            
            break;
        }
        
        else
        {
            std::cout << "Greska! Pokusajte ponovo." << std::endl;
        }
    }
    
    return 0;
}
