#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string.h>

using namespace std;
const double accuracy = 0.001; // for interpolation

const double toRad = 0.017453293;
double const c = 3e+8;           // speed of light
double const u = 931.43432;      // mass unit
double const aum = 931494.32; // MeV/c^2
string srimFiles[2];
struct Data
{
    int mass[4];
    int charge[4];
    string ionName[4];
    string fmass;
  	double eBeam;
    double thickness;
    double distanceFromWindow;
    double labAngle;
    
    double distanceToOR;
     
};
void trancateFile(string input, string output)
{

 ifstream from(input.c_str());
 if (!from) cout<<"Cannot open input file";

 ofstream to(output.c_str());
 if (!to) cout<<"Cannot open output file";


 char ch;
 char buffer[1000];

 while(from.get(ch)){
    if(ch=='S'){
     
    from.get(ch);
    if(ch=='t'){
      from.get(ch);
      if(ch=='r'){

          from.get(ch);
          if(ch=='a'){
             from.getline(buffer,1000);
              from.getline(buffer,1000);
              
          while(from.get(ch)) {
                 if(ch == '-') {
                if(from.peek() == '-') break;
                }   
                while(from.peek()==','){
                  from.ignore(1,',');
                    from.putback('.');
                }
                  to.put(ch);  
                  }
          }
          else{
            continue;
          }
      }
      else{
        continue;
      }
    }
      else{
        continue;
      }
    }
   }
   





cout <<"The file "<<"\""<< output <<"\""<< " was created" << endl;

from.close();
to.close();

}
void toTwocolumn(string fromName, string toName)
{
  ifstream from;
  ofstream to;

  char buf[100];
  char unit[3];
  double energy,eloss;
  int iter(0);
 
    
  from.open(fromName.c_str());
  if(!from)cout<<"Cannot open input file"; 
     
  to.open(toName.c_str());
  if(!to) cout<<"Cannot open input file";
  
  to.setf(ios::fixed);
  to.setf(ios::showpoint);
  to.precision(3);


  while(1) {
      from >> energy >> unit >> eloss; from.getline(buf,100);
        if (strcmp(unit,"keV") == 0) 
           energy /=1000;
        to << setw(8) << energy << setw(8) << eloss << endl;
        iter++;
          if (!from.good()) break;
  }
  // cout << iter << " good lines were read" << endl;
    
  from.close();
  to.close();
  cout <<"The file "<<"\""<< toName <<"\""<< " was created" << endl;
}
double interpolation(int k, double energyCur, vector<double> e, vector<double> stopping_power)
{
  double var  =  energyCur;
  double e1 =  e.at(k-2);
  double e2 =  e.at(k-1);
  double e3 =  e.at(k);
  double dr1 = stopping_power.at(k-2);
  double dr2 = stopping_power.at(k-1);
  double dr3 = stopping_power.at(k);
  
  double a,b,c;

 /*
 Lagrange's classical forumula of polynomial interpolation.
 I used 3 points interpolation, while it can be used more.
*/

 a = (var-e2)*(var-e3)/((e1-e2)*(e1-e3))*dr1;
 b = (var-e1)*(var-e3)/((e2-e1)*(e2-e3))*dr2;
 c = (var-e1)*(var-e2)/((e3-e1)*(e3-e2))*dr3;
 
 return (a+b+c);

}
void simpsonInteg(string fromName, string toName, int &number)
{
  double eCurrent;
  double eStart;
  double de(0),dedr(0);

  int index(0);
  int j(0);
  double range(0);
  double emin;
  double step;

  vector<double> e;  
  vector<double> stopping_power;

  double t1,t2;
  
  ifstream from;
  ofstream to;
    
  from.open(fromName.c_str());
  if(!from) cout<<"Cannot open input file"<<endl; 
     
  to.open(toName.c_str());
  if(!to) cout<<"Cannot open output file"<<endl;  
    
  to.setf(ios::fixed);
  to.setf(ios::showpoint);
  to.precision(4);
  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(5);

 
  
  number = 0;

   while(1) {
     from >> t1 >> t2;
     e.push_back(t1);
     stopping_power.push_back(t2);
     j++;
       if (!from.good()) break;
   }

  index = j - 2;
  

  eCurrent = e.at(index);
  eStart   = e.at(index);

  emin = e.at(1);
  step = accuracy*e.at(index)/stopping_power.at(index);
  step = accuracy;
 

  while (eCurrent > emin) 
  {
     if(eCurrent <= e.at(index-1))
       index--;

      dedr = interpolation(index,eCurrent,e,stopping_power);
      step = accuracy*eCurrent/stopping_power.at(index);   

      de = dedr*step;
      range+=step;
      eCurrent-=de;
     
      to << setw(10) << range << setw(18) << dedr
         << setw(22)<< (eStart -eCurrent)*1000 << endl;   

     number++;
      
  }

 // to << setw(10) << range << setw(18) << dedr
  //   << setw(22) << eStart*1000 << endl; 

  number+=1;

 from.close();
 to.close();
 cout <<"The file "<<"\""<< toName <<"\""<< " was created" << endl;
}
double timeFlight(double e1, double e2, double m, double x)
{
 
double vi = c*sqrt(2*e1/(m*u));
double vf = c*sqrt(2*e2/(m*u));
if(vf<1){

  return x/vi;
}
if(vf!=vi){
double acc=(vi*vi-vf*vf)/2/x;//acceleration
return (vi-vf)/acc;       
}
else{
 return x/vi;
}

}
double qReaction(int *im, int *z, double *dm, string fp)
{
   int iter(0);
   int charge,mass;
   double defmass,errdef;
   
   double ddef[4];
 
   
   ifstream from(fp);
   if(!from) cout<<"Cannot open input file ---- "<<endl; 
     
  
   while (iter < 4 && !from.eof())
   {
     from >> charge >> mass >> defmass >> errdef;
     
      for (int i = 0; i < 4; i++)
      {
        if (charge == z[i] && mass == im[i])
	{
	      ddef[i] = defmass;
	      iter++;	      
        }
      }
    }
    
      for (int i = 0; i < 4; i++)
        dm[i]= ddef[i]/aum + (double)im[i];  
    
    
    from.close();
      
  return (ddef[0] + ddef[1] - ddef[2] - ddef[3])*0.001;
}
//________

void kinTwoBody(double *m, double beamEn, double eex, double q, double phirad, double &e3_first, double &e3_second,double &e4_first,double &e4_second, double &e3cm, double &int1, double &int2, bool &key)
{

 double phi_max;
 double a,b,c,d;

 double q_tot;
 double e_tot;
 double temp;
 
 q_tot = q - eex;
 e_tot = beamEn + q_tot;	 

 a = m[0]*m[3]*(beamEn/e_tot)/(m[0]+m[1])/(m[2]+m[3]);
 b = m[0]*m[2]*(beamEn/e_tot)/(m[0]+m[1])/(m[2]+m[3]);
 c = m[1]*m[2]/(m[0]+m[1])/(m[2]+m[3])*(1 + m[0]*q_tot/m[1]/e_tot);
 d = m[1]*m[3]/(m[0]+m[1])/(m[2]+m[3])*(1 + m[0]*q_tot/m[1]/e_tot);

 //cout.precision(3);

 key = false;

 temp = sqrt(d/b - pow(sin(phirad),2));

 if (b <= d) 
 { 
     e3_first = e_tot*b*pow(cos(phirad)+temp,2);
     e3_second = 0.0;
    e4_first=e_tot-e3_first;
    e4_second=0.0;
 } 
 else 
 {
        key = true;
        phi_max = asin(sqrt(d/b));

        if(phi_max >= phirad)
	{
	   e3_first  = e_tot*b*pow(cos(phirad)+temp,2);
           e3_second = e_tot*b*pow(cos(phirad)-temp,2);
       e4_first=e_tot-e3_first;
    e4_second=e_tot-e3_second;     
           //error("Two branches of energy is possible","");
	}
	else
	{
	   e3_first = 0.0;
           e3_second = 0.0;
           cout<<"Lab angle of light particle      : "<< phirad/toRad << endl;
           cout<<"The maximum angle from kinematics: "<< phi_max/toRad << endl;
           //error("Two branches of energy is possible. Limit Angle was reached","");
         
        }
 }

   e3cm = d*e_tot;
   
   int1 = 0;
   int2 = 0;

    if (e3_first != 0)
      int1 = sqrt(a*c)*temp/(e3_first/e_tot);

    if (e3_second != 0)
      int2 = sqrt(a*c)*temp/(e3_second/e_tot);
    
}
int readElements(string name, vector<double> &par1, vector<double> &par2, vector<double> &par3)
{
  int number = 0;
  ifstream from;
  
  double t1,t2,t3;
  
  from.open(name);
         
  
   
	while(1) {	
         if (!from.good()) break;    
           from >> t1 >> t2 >> t3;
           par1.push_back(t1);
           par2.push_back(t2);
           par3.push_back(t3);
        }

  number = par1.size();

  from.close();
  return number-1;
}
double elosSearch(double range, double current, vector<double> &ar1, vector<double> &ar2, vector<double> &ar3, int size, double &dedr)
{
     double emax = ar3.at(size-1);
     double e_more(0), e_less(0);
     int iter(0);
     double x0,x1;
     double elow, elos;
     
      if(current > emax) {
         cout <<"ECurrent " << setw(16) << current*0.001 << " MeV" << endl;
         cout <<"EMax [from SRIM] " << emax*0.001 << " MeV" << endl;
          cout <<"You need to increase your  EMax [from SRIM] "<< endl;
       return 0;
     }

      while(1)
      {
        e_less = emax - ar3.at(iter);
        iter++;
          if(e_less < current) break;
        e_more = e_less; 
      }

     if(iter == 1)
       x0 = ar1.at(iter-1)*(emax-current)/ar3.at(iter-1);
     else
       x0 = ar1.at(iter-2)+(ar1.at(iter-1)-ar1.at(iter-2))*(1-(current-e_less)/(e_more-e_less)); 

      
     x1 = x0 + range;
     
     if(x1 > ar1.at(size-1)) {
        return current*0.001;
     }

     while(x1 > ar1.at(iter-1))  
      iter++;

      if(iter - 2 < 0)
      {
         dedr = ar2.at(iter-1); 
         return 0;    
      }
   
     elow = emax-(ar3.at(iter-2)+(x1-ar1.at(iter-2))/(ar1.at(iter-1)-ar1.at(iter-2))*(ar3.at(iter-1)-ar3.at(iter-2)));
     elos = current - elow;        
     dedr = ar2.at(iter-1);     
 
  return elos/1000;           
}
void gobble(ifstream& f, char ch)
{
 while(f.peek() == ch)
  while(f.get() != '\n')
    {
    }
}
void readParam(std::ifstream& f,char *iName, Data &d,  string fSrim[])
{
 
  char buffer[100];
  int size=2;
   f.open(iName);
   
   if(!f) (std::cout<<"Cannot open input file") ;
     
   else {
  
   gobble(f,'#');
		   

    for(int j=0; j < 4; j++) {
      f >> d.mass[j] >> d.charge[j] >> d.ionName[j];
      f.getline(buffer,100);
   }

   f >> d.eBeam;
   f.getline(buffer,100);
  
   f >> d.thickness;
   f.getline(buffer,100);
   f >> d.distanceFromWindow;
   f.getline(buffer,100);
   f >> d.distanceToOR;
   f.getline(buffer,100);
   f >> d.labAngle;
   f.getline(buffer,100);
   gobble(f,'#');
   
  
   gobble(f,'#');
   
   f >> fSrim[size-2];
   f.getline(buffer,100);
   f >> fSrim[size-1];
   f.getline(buffer,100);

   gobble(f,'#');




   f >> d.fmass;
   f.getline(buffer,100);

   f.close();

   std::cout <<"Information from file "<<"\"" <<iName <<"\""<<" was read"<< endl;

 }       
}
int main(){

  double ecm, ecm_recoil, exener;         	 
  double e_proj, eBeam_current; 
  bool key;
  double qOfReaction;
  double doubleMass[4];

  double phiRad,ry,rx;
  string fmass;
  vector<double> rangeHeavy,rangeLight;
  vector<double> dedrHeavy,dedrLight;
  vector<double> deHeavy,deLight;
   int heavyNum, lightNum;

  int heavySize, lightSize;

   string srimTrancated[2];
   string srimTwoCol[2];
   string srimThreeCol[2];
   
  double stop_pow_heavy, stop_pow_light;
  double elos_beam;
  double heavyMaxEnergy,lightMaxEnergy,heavyMinEnergy;
  double range(0),step,range_left(0), total_range(0);  
  ofstream TOF;
  

  double edet,edet_old(0),coef,zcan,kdedx,omega(0);  
  double e1,e2,de1,de2,elos1,elos2,eEx(0);
  double ks1,ks2,ks,e3_first,e3_second,e4_first,e4_second;
  double time1, time2, dtime;
  ifstream paramFile;

  Data expData;
  char* input="N.dat";


  readParam(paramFile,input,expData,srimFiles);



  bool choice=true;     //true if elastic pro or direct reaction      


  
 
  int count(0);
  
  
  string foutput="tof"+expData.ionName[2]+".dat";
  TOF.open(foutput);   
  
   TOF << "#---------------------  #\n";
   TOF << "#| Energy in  |  Time  |#\n";
   TOF << "#|  Detector  |        |#\n";
   TOF << "#|   'MeV'    |  'ns'  |#\n";
   TOF << "#-----------------------#\n";

    srimTrancated[0]   = expData.ionName[0]+".Tred";
   srimTrancated[1]   = expData.ionName[2]+".Tred";
   srimTwoCol[0]      = expData.ionName[0]+".2col";
   srimTwoCol[1]      = expData.ionName[2]+".2col";
   srimThreeCol[0]    = expData.ionName[0]+".3col";
   srimThreeCol[1]    = expData.ionName[2]+".3col";

  phiRad = toRad*expData.labAngle;
  double det_angle = phiRad;

  trancateFile(srimFiles[0], srimTrancated[0]);
  trancateFile(srimFiles[1], srimTrancated[1]);
 
   // Creating files from original srim files
   toTwocolumn(srimTrancated[0],srimTwoCol[0]);     
   toTwocolumn(srimTrancated[1],srimTwoCol[1]);

   // Making interpolation of srim files
   simpsonInteg(srimTwoCol[0],srimThreeCol[0],heavyNum);  
   simpsonInteg(srimTwoCol[1],srimThreeCol[1],lightNum);


  ry = expData.distanceToOR*sin(det_angle);
  rx = expData.distanceToOR*cos(det_angle) + expData.distanceFromWindow;
  total_range = expData.thickness*(expData.distanceFromWindow + expData.distanceToOR);
     qOfReaction=qReaction(expData.mass,expData.charge,doubleMass,expData.fmass);


  heavySize = readElements(srimThreeCol[0],rangeHeavy,dedrHeavy,deHeavy);
  heavyMaxEnergy = deHeavy.at(heavySize-1); 
  heavyMinEnergy = deHeavy.at(0); 

   lightSize = readElements(srimThreeCol[1],rangeLight,dedrLight,deLight);
  lightMaxEnergy = deLight.at(lightSize-1); 
  range = 0;
  step = total_range/1000; 
while (1) 
 {
    if(range > total_range) {
       // cout << range <<"  " << total_range << endl;
      cout<<"Correct thickness of target, it is too small"<<endl;
        cout <<"The file "<<"\""<<  " was created" << endl;
      return 0;
   }
  
   // cout << range <<"  " << total_range << endl;
  
     phiRad = atan(ry/(rx - range/expData.thickness));
    
    
   // cout << ry << " " << rx << " " << phiRad << " " << total_range << endl;

    count++;

   if(count == 1 ) 
   {
     elos_beam = elosSearch(range,expData.eBeam*1000,rangeHeavy,dedrHeavy,deHeavy,heavySize,stop_pow_heavy);
     eBeam_current = expData.eBeam;
     elos_beam = 0.0;
   }
   else
   {
     elos_beam = elosSearch(range,expData.eBeam*1000,rangeHeavy,dedrHeavy,deHeavy,heavySize,stop_pow_heavy);
     eBeam_current = expData.eBeam - elos_beam;
   }   
   
   //cout << eBeam_current << " " << heavyMinEnergy/1000 <<endl;
   
   if(eBeam_current < heavyMinEnergy/1000) {
      cout <<"The file "<<"\""<<  " was created" << endl;
      return 0;
   }
 

    kinTwoBody(doubleMass,eBeam_current,eEx,qOfReaction,phiRad,e3_first,e3_second,e4_first,e4_second,ecm_recoil,ks1,ks2,key);
  //  cout << e3_first << " " << e3_second <<endl;  


 // Conditions to select elastic or inelastic treatment
 
   if(key && ((e3_first || e3_second) < 0.0001)){
      cout << "Two branches are possible. Limit angle achieved" << endl;    
      return 0;
     }

    if(choice) {
     e_proj = e3_first;
     ks = ks1;
   //  cout << choice <<endl;
    }
    else {
     e_proj = e3_second;
     ks = ks2;
     //  cout << choice <<endl;
     
    }   

  
//   kinTwoBody(m,eBeam_current,d.eEx,q,phiRad,e_proj,ecm_recoil,ks1);
    
   if (e_proj < 0.001) {
        cout <<"Warning: The projectile energy is " << e_proj <<endl;
        return 0;
     }
   
   coef = eBeam_current/e_proj;

     range_left = (rx - range/expData.thickness)*expData.thickness;
  
  
   e1 = e_proj;

   elos1 = elosSearch(range_left,e1*1000,rangeLight,dedrLight,deLight,lightSize,stop_pow_light);
   
   de1 = e1 - elos1;
 
  
 

   
    time1 = timeFlight(expData.eBeam,eBeam_current,doubleMass[0],range/expData.thickness*0.01);
    time2 = timeFlight(e1,abs(de1),doubleMass[2],range_left/expData.thickness*0.01);          ///!!!
    dtime = time1 + time2; 
  edet = de1;

  
   //if(elos1 > e1) break; eto bilo tak iznachal'no mnoi napisano

  range+=step; 

   
  
   cout << setw(10) << range/expData.thickness  <<setw(14) <<eBeam_current<< setw(14)<<e1 << setw(14)<< de1 << setw(10)<< dtime*1e9  <<  endl;
  
  TOF<< setw(11)<< de1 <<setw(11)<< dtime*1e9  <<  endl;
 
}
TOF.close();

}