//Author: Gabriele Inghirami (g.inghirami@gsi.de)
//License: PUBLIC DOMAIN 
//References in preprints on arxiv.org: 1602.02223, 1305.5806, 1609.03042

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

const double h=0.197326;
const double sqrh=0.444214;//it is sqrt(hbar*c)=sqrt(0.197326)
const double s=0.0058;//electrical conductivity
const double schir=0.0015;//chiral conductivity
const double collision_energy=200;//sqrt(s) of the collision, in GeV/nucleon
const double Qbig=0.30282212;//it is sqrt(4*Pi*alpha)
const double R=6.38;//Au
//const double R=6.62;//Pb
double B; //impact parameter in fm
const double mass_number=197;//Au
//const double mass_number=208;//Pb
const double Zcharge=79;//number of protons of the colliding ion Au
//const double Zcharge=82;//number of protons of the colliding ion Pb
double g,v; //gamma Lorentz factor and ion velocity
const int coordinates=2; //1=Minkowski, 2=Milne
char* outputfile;
const int printstdout=0;//1=print results of integrations to std output, 0=don't print
const int theta_active=1; //we suppress the contribution for z-zlimit>v*t

const double gsl_relerr=1.e-3; //it controls the relative error in integration
const double gsl_abserr=1.e-6; //it controls the absolute error in integration
const double coarse_error_factor1=10;//if the first integration fails, we retry using a tolerance larger by this factor
const double coarse_error_factor2=100;//if the second integration fails, we retry using a tolerance larger by this factor
const double intervals_limit=100000; //maximum number of subintervals

/* data grid */
const int xcells=141;
const int ycells=141;
int zcells;
const double xside=14;
const double yside=14;
double* zarr;
const double t0=0.4;

int failed; //flag to mark the integrations with errors

/* projection contains informations about the nucleus and the field component */
typedef struct f_projection{int sign; double (*trigo)(const double, const double, const double, const double);} projections;

/* we define a general type for parameters, valid for all cases, but each time we´ll fill it with the needed values only */
typedef struct f_params{ double t; double x; double y; double z; double xi; double yi; int nucleus; int component; int Btype; projections projection;
	                     double (*dist)(const double, const double); double (*zpos)(const double, const double, const double);} parameters;

void printhelp(int);

double H0(parameters *,int, int, int); //second parameter is for nucleus 1 or 2, third parameter for x and y components, fourth parameter for the B field type (0=all, 1=classic, 2=chiral)
double H0_integrand(double, void *);

double H1_internal(parameters *);
double H1_internal_integrand(double, void *);


double Hpoint(parameters *);

double dist1(const double, const double);
double dist2(const double, const double);

double zpos1(const double, const double, const double);
double zpos2(const double, const double, const double;);

double Cpsi(const double, const double, const double, const double);
double Spsi(const double, const double, const double, const double);
double Noangle(const double, const double, const double, const double);

int main(int argc,char *argv[])
{
    parameters input_params;
    double result, uncertainty;
    double H1x, H2x, H1y, H2y, H1Cx, H2Cx, H1Cy, H2Cy, H1Cz, H2Cz;
    double dx, xcoord, dy, ycoord, zcoord,dz,xmin,ymin;
    int x,y,z,xxcells,yycells;
    double zmin, zmax;
    int npoints;
    double coord_conversion_factor, z_conversion_factor;
    int Borigin; //0=both chiral and classic, 1=classic only, 2=chiral only
    FILE* fout;
    int H1x_on, H2x_on, H1y_on, H2y_on, H1Cx_on, H2Cx_on, H1Cy_on, H2Cy_on, H1Cz_on, H2Cz_on;
    int qq; //just a counter
    const char resok[]="OK";
    const char resno[]="NO";
    char* resmesg[10];
    int failcode[10];
    int printfile;
    int kind_of_run;//1=sequence, 2=range, 3=point
    int protons, neutrons;
    double average_nucleon_mass,yb;


    if(argc < 6)
    {
       print_syntax_and_exit();
    }

    if(strncmp(argv[2],"seq",3)==0)
    {
      if(argc < 6)
      {
       print_syntax_and_exit();
      }
      else
      {
       kind_of_run=1;
       printfile=1;
       xxcells=xcells;
       yycells=ycells;
       xmin=-xside;
       ymin=-yside;
      }
    }
    else if(strncmp(argv[2],"range",5)==0)
    {
      if(argc != 8)
      {
       print_syntax_and_exit();
      }
      else
      {
       kind_of_run=2;
       printfile=1;
       xxcells=xcells;
       yycells=ycells;
       xmin=-xside;
       ymin=-yside;
      }
    }
    else if(strncmp(argv[2],"point",5)==0)
    {
      if(argc != 7)
      {
       print_syntax_and_exit();
      }
      else
      {
       kind_of_run=3;
       printfile=0;
       xxcells=1;
       yycells=1;
       xmin=-atof(argv[4]);
       ymin=-atof(argv[5]);
      }
    }
    else
    {
       print_syntax_and_exit();
    }
       
    Borigin=atoi(argv[3]);

    if((Borigin < 0) || (Borigin >2)) {
      printf("Erron in parameter Borigin, which can be only 0,1,2\n");
      exit(2);
    }

    if(kind_of_run==1)
    {
      outputfile=argv[4];
      zcells=argc-4;

      zarr=malloc(zcells*sizeof(double));
      if(zarr == NULL) {
        printf("Error in allocating the zarr array. Leaving.\n");
        exit(2);
      }

      for(qq=0;qq<zcells;qq++)
      {
        zarr[qq]=atof(argv[4+qq]);
      } 
    }
    else if(kind_of_run==2)
    {
      outputfile=argv[4];
      zmin=atof(argv[5]);
      zmax=atof(argv[6]);
      if(zmin>zmax) {
        printf("Error, zmin cannot be bigger than zmax!\n");
        exit(2);
      }
      zcells=atoi(argv[7]);
      if(zcells<1) {
        printf("Error, the number of z cells must be > 0!\n");
        exit(2);
      }
      dz=(zmax-zmin)/zcells;

      zarr=malloc(zcells*sizeof(double));
      if(zarr == NULL) {
        printf("Error in allocating the zarr array. Leaving.\n");
        exit(2);
      }

      for(qq=0;qq<zcells;qq++)
      {
        zarr[qq]=zmin+(0.5+qq)*dz;
      } 
    }

    B=atof(argv[1]);
    if(B<0) 
    { 
       printf("Error, the impact parameter B must be > 0!!!\n");
       exit(1);
    }

    dx=(xside*2)/xxcells;
    dy=(yside*2)/yycells;

    protons=Zcharge;
    neutrons=(mass_number-protons);
    average_nucleon_mass=(protons*0.9382721+neutrons*0.9395654)/mass_number;

    yb=log(collision_energy/average_nucleon_mass);

    g=collision_energy/(2*average_nucleon_mass);

    v=sqrt(1. - 1/(g*g));

    if(!((coordinates == 1) || (coordinates == 2)))
    {
     printf("coordinates unproperly set; it can be only 1 or 2\n");
     exit(3);
    }
    
    if(kind_of_run == 3) {
       dx=0;
       dy=0;
       zcells=1;
       printfile=0;
       zarr=malloc(zcells*sizeof(double));
       if(zarr == NULL) {
         printf("Error in allocating the zarr array. Leaving.\n");
         exit(2);
       }
       zarr[0]=atof(argv[6]);
    }
                
    if(printfile==1) {
      fout=fopen(outputfile, "w");
      if(fout == NULL) {
        printf("Error in opening the output file %s. I quit.\n",outputfile);
        exit(1);
      }
      fprintf(fout,"#B field for: b_impact: %lf, tau: %lf, coord: %d, Sqrt(sNN): %lf, nucleons: %d, protons: %d, nucl. radius: %lf, sigma: %lf, sigma_chiral: %lf\n", B, t0, coordinates, collision_energy, (int) mass_number, (int) Zcharge, R, s, schir);
      fprintf(fout,"#1: x  2: y  3: z or eta  4: H1x  5: msg+failcode  6: H2x  7: msg+failcode   8: H1Cx  9: msg+failcode   10: H2Cx  11: msg+failcode   12: H1y  13: msg+failcode  14: H2y  15: msg+failcode  16: H1Cy  17: msg+failcode  18: H2Cy  19: msg+failcode  20: H1Cz  21: msg+failcode   22: H2Cz  23: msg+failcode\n"); 
    } 
    if((printstdout==1) || (kind_of_run==3)) printf("#1: x  2: y  3: z or eta  4: H1x  5: msg+failcode  6: H2x  7: msg+failcode   8: H1Cx  9: msg+failcode   10: H2Cx  11: msg+failcode   12: H1y  13: msg+failcode  14: H2y  15: msg+failcode  16: H1Cy  17: msg+failcode  18: H2Cy  19: msg+failcode  20: H1Cz  21: msg+failcode   22: H2Cz  23: msg+failcode\n");
 
    for(z=0;z<zcells;z++)
       {
          zcoord=zarr[z];

          if(coordinates == 1)
          {
            input_params.t=t0;
            input_params.z=zcoord;
            coord_conversion_factor=1;
            z_conversion_factor=1;
          }
          else
          {
            input_params.t=t0*cosh(zcoord);
            input_params.z=t0*sinh(zcoord);
            coord_conversion_factor=1./cosh(zcoord);
            z_conversion_factor=1./t0;
          }

          for(y=0;y<yycells;y++)
          {
           input_params.y=ymin+(y+0.5)*dy;
           for(x=0;x<xxcells;x++)
              {
                input_params.x=xmin+(x+0.5)*dx; 
                H1x=0;
                H2x=0;
                H1y=0;
                H2y=0;
                H1Cx=0;
                H2Cx=0;
                H1Cy=0;
                H2Cy=0;
                H1Cz=0;
                H2Cz=0;
                H1x_on=1;
                H2x_on=1;
                H1y_on=1;
                H2y_on=1;
                H1Cx_on=1;
                H2Cx_on=1;
                H1Cy_on=1;
                H2Cy_on=1;
                H1Cz_on=1;
                H2Cz_on=1;
                //By default everything is active, now we suppress some computations
                if(Borigin == 1)
                {
                  H1Cx_on=0;
                  H2Cx_on=0;
                  H1Cy_on=0;
                  H2Cy_on=0;
                  H1Cz_on=0;
                  H2Cz_on=0;
                }
                if(Borigin == 2)
                {
                  H1x_on=0;
                  H2x_on=0;
                  H1y_on=0;
                  H2y_on=0;
                }
                if(input_params.y == 0) {  
                  H1x_on=0;
                  H2x_on=0;
                  H1Cy_on=0;
                  H2Cy_on=0;
                  }
                if(2*input_params.x == -B) {
                  H1y_on=0;
                  H1Cx_on=0;
                  }
                if(2*input_params.x == B)  {
                  H2y_on=0;
                  H2Cx_on=0;
                  }

                if(theta_active == 1) {
                  if(zpos1(input_params.t,input_params.z,v)<=0)
                    {
                     H1x_on=0;
                     H1y_on=0;
                     H1Cx_on=0;
                     H1Cy_on=0;
                     H1Cz_on=0;
                     resmesg[0]=resok;
                     resmesg[2]=resok;
                     resmesg[4]=resok;
                     resmesg[6]=resok;
                     resmesg[8]=resok;
                  }
                  if(zpos2(input_params.t,input_params.z,v)<=0)
                    {
                     H2x_on=0;
                     H2y_on=0;
                     H2Cx_on=0;
                     H2Cy_on=0;
                     H2Cz_on=0;
                     resmesg[1]=resok;
                     resmesg[3]=resok;
                     resmesg[5]=resok;
                     resmesg[7]=resok;
                     resmesg[9]=resok;
                  }
                }
                     
                if(H1x_on==1) 
                {
                   failed=0;
                   failcode[0]=0;
                   H1x=H0(&input_params,1,1,1)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[0]=resok;
                   } else
                   {
                    resmesg[0]=resno;
                   } 
                   failcode[0]=failed;
                } 
                else
                {
                   resmesg[0]=resok;
                   failcode[0]=0;
                } 

                if(H2x_on==1)
                {
                   failed=0;
                   failcode[1]=0;
                   H2x=H0(&input_params,2,1,1)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[1]=resok;
                   } else
                   {
                    resmesg[1]=resno;
                   } 
                   failcode[1]=failed;
                } 
                else
                {
                   resmesg[1]=resok;
                   failcode[1]=0;
                } 

                if(H1Cx_on==1)
                {
                   failed=0;
                   failcode[2]=0;
                   H1Cx=H0(&input_params,1,1,2)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[2]=resok;
                   } else
                   {
                    resmesg[2]=resno;
                   } 
                   failcode[2]=failed;
                } 
                else
                {
                   resmesg[2]=resok;
                   failcode[2]=0;
                } 

                if(H2Cx_on==1)
                {
                   failed=0;
                   failcode[3]=0;
                   H2Cx=H0(&input_params,2,1,2)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[3]=resok;
                   } else
                   {
                    resmesg[3]=resno;
                   } 
                   failcode[3]=failed;
                } 
                else
                {
                   resmesg[3]=resok;
                   failcode[3]=0;
                } 

                if(H1y_on==1)
                {
                   failed=0;
                   failcode[4]=0;
                   H1y=H0(&input_params,1,2,1)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[4]=resok;
                   } else
                   {
                    resmesg[4]=resno;
                   } 
                   failcode[4]=failed;
                } 
                else
                {
                   resmesg[4]=resok;
                   failcode[4]=0;
                } 

                if(H2y_on==1)
                {
                   failed=0;
                   failcode[5]=0;
                   H2y=H0(&input_params,2,2,1)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[5]=resok;
                   } else
                   {
                    resmesg[5]=resno;
                   } 
                   failcode[5]=failed;
                } 
                else
                {
                   resmesg[5]=resok;
                   failcode[5]=0;
                } 

                if(H1Cy_on==1)
                {
                   failed=0;
                   failcode[6]=0;
                   H1Cy=H0(&input_params,1,2,2)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[6]=resok;
                   } else
                   {
                    resmesg[6]=resno;
                   } 
                   failcode[6]=failed;
                } 
                else
                {
                   resmesg[6]=resok;
                   failcode[6]=0;
                } 

                if(H2Cy_on==1)
                {
                   failed=0;
                   failcode[7]=0;
                   H2Cy=H0(&input_params,2,2,2)*coord_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[7]=resok;
                   } else
                   {
                    resmesg[7]=resno;
                   } 
                   failcode[7]=failed;
                } 
                else
                {
                   resmesg[7]=resok;
                   failcode[7]=0;
                } 

                if(H1Cz_on==1)
                {
                   failed=0;
                   failcode[8]=0;
                   H1Cz=H0(&input_params,1,3,2)*z_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[8]=resok;
                   } else
                   {
                    resmesg[8]=resno;
                   } 
                   failcode[8]=failed;
                } 
                else
                {
                   resmesg[8]=resok;
                   failcode[8]=0;
                } 

                if(H2Cz_on==1)
                {
                   failed=0;
                   failcode[9]=0;
                   H2Cz=H0(&input_params,2,3,2)*z_conversion_factor;
                   if(failed == 0) 
                   {
                    resmesg[9]=resok;
                   } else
                   {
                    resmesg[9]=resno;
                   } 
                   failcode[9]=failed;
                } 
                else
                {
                   resmesg[9]=resok;
                   failcode[9]=0;
                } 

                if(printfile==1) fprintf(fout,"%7.5f %7.5f %7.5f %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d\n", input_params.x, input_params.y,  zcoord, H1x, resmesg[0],failcode[0], H2x, resmesg[1],failcode[1], H1Cx, resmesg[2],failcode[2], H2Cx, resmesg[3],failcode[3], H1y, resmesg[4],failcode[4], H2y, resmesg[5],failcode[5], H1Cy, resmesg[6],failcode[6], H2Cy, resmesg[7],failcode[7], H1Cz, resmesg[8],failcode[8], H2Cz, resmesg[9],failcode[9]);  
                if((printstdout==1) || (kind_of_run==3)) printf("%7.5f %7.5f %7.5f %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d %g %s%02d\n", input_params.x, input_params.y,  zcoord, H1x, resmesg[0],failcode[0], H2x, resmesg[1],failcode[1], H1Cx, resmesg[2],failcode[2], H2Cx, resmesg[3],failcode[3], H1y, resmesg[4],failcode[4], H2y, resmesg[5],failcode[5], H1Cy, resmesg[6],failcode[6], H2Cy, resmesg[7],failcode[7], H1Cz, resmesg[8],failcode[8], H2Cz, resmesg[9],failcode[9]);  
              }
           if(printfile==1) fprintf(fout,"\n");
          }
          if(printfile==1) fprintf(fout,"\n");
       } 
    if(printfile==1)
    {
     fclose(fout);
    }
  
    return 0;
}

 
double Hpoint(parameters *Hparams)
{
    double uncertainty;
    int error;
    double t = (Hparams->t);
    double x = (Hparams->x);
    double y = (Hparams->y);
    double z = (Hparams->z);
    double xi= (Hparams->xi);
    double yi= (Hparams->yi);
    double zpos_point;
    double xt2,xt,D,sqrD,A;
   
    zpos_point=Hparams->zpos(t,z,v);
    
    xt2=(x-xi)*(x-xi)+(y-yi)*(y-yi);
    xt=sqrt(xt2);
    D=g*g*zpos_point*zpos_point+xt2;
    sqrD=sqrt(D);
    A=(s*v*g/2)*(g*zpos_point-sqrD)/h;

    if(Hparams->Btype == 1) {
      return (xt/pow(D,1.5))*(1+s*v*g*sqrD/(2*h))*exp(A);
    } else {
      if(Hparams->component==3) {//chiral case, Bz
        return (g*g*zpos_point*zpos_point*(1+s*v*g*sqrD/(2*h))+D*(1-s*v*g*sqrD/(2*h)))*exp(A)/pow(D,1.5);
      } else {//chiral case, Bx or By
        return (xt/pow(D,1.5))*(g*zpos_point+A*sqrD)*exp(A);
      }
    }

}

double H1_internal(parameters *Hparams)
{
    double result, uncertainty;
    double yl, yr; //integration limits
    
    double xi=Hparams->xi;
    int err_gsl;
    
    if(Hparams->nucleus==1)
    {
		yr=sqrt(R*R - (xi + B/2.)*(xi + B/2.));
    }
    else //we already cheched that nucleus can be only 1 or 2
    {
		yr=sqrt(R*R - (xi - B/2.)*(xi - B/2.));
    }
	yl=-yr;

    gsl_set_error_handler_off();
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (intervals_limit);
  
    gsl_function F;
    F.function = &H1_internal_integrand;
    F.params = Hparams;
 
    err_gsl=gsl_integration_qags(&F, yl, yr, gsl_abserr, gsl_relerr, intervals_limit, w, &result, &uncertainty);
 
    if(err_gsl != 0) {
      failed=1;
      err_gsl=gsl_integration_qags(&F, yl, yr, gsl_abserr*coarse_error_factor1, gsl_relerr*coarse_error_factor1, intervals_limit, w, &result, &uncertainty);
      if(err_gsl != 0) {
        failed=2;
        err_gsl=gsl_integration_qags(&F, yl, yr, gsl_abserr*coarse_error_factor2, gsl_relerr*coarse_error_factor2, intervals_limit, w, &result, &uncertainty);
        if(err_gsl !=0) failed=3;
      }
    }
	
    gsl_integration_workspace_free (w);


    return result;	
}
 
double H0(parameters *Hparams, int nucleus, int component, int Btype)
{
    double result, uncertainty;
    double xl, xr; /* left and right integration boundaries along x */
    double rho;
    int integration_error;
    
    rho=Zcharge/(4./3.*M_PI*pow(R,3.));
        
    Hparams->nucleus=nucleus;
    Hparams->component=component;
    Hparams->Btype=Btype;


    if(nucleus == 1)
    {
		xl=-R-B/2.;
		xr=R-B/2.;
		Hparams->dist=&dist1;
		Hparams->zpos=&zpos1;
                if(Btype == 1) {
		  switch(component)
		  {
			case 1: Hparams->projection.trigo=&Spsi; Hparams->projection.sign=+1; break;
			case 2: Hparams->projection.trigo=&Cpsi; Hparams->projection.sign=-1; break;
	                default: printf("Component is unclear, leaving...\n"); exit(2); break;
		  }
                }
                else {
		  switch(component)
		  {
			case 1: Hparams->projection.trigo=&Cpsi; Hparams->projection.sign=+1; break;
			case 2: Hparams->projection.trigo=&Spsi; Hparams->projection.sign=+1; break;
                        //we have the case 3 only with chiral B field
                        case 3: Hparams->projection.trigo=&Noangle; Hparams->projection.sign=-1; break;
	                default: printf("Component is unclear, leaving...\n"); exit(2); break;
		  }
                }
	}
	else if(nucleus == 2)
	{
		xl=-R+B/2.;
		xr=R+B/2.;
		Hparams->dist=&dist2;
		Hparams->zpos=&zpos2;
                if(Btype == 1) {
		  switch(component)
		  {
			case 1: Hparams->projection.trigo=&Spsi; Hparams->projection.sign=-1; break;
			case 2: Hparams->projection.trigo=&Cpsi; Hparams->projection.sign=+1; break;
	                default: printf("Component is unclear, leaving...\n"); exit(2); break;
		  }
                }
                else{
		  switch(component)
		  {
			case 1: Hparams->projection.trigo=&Cpsi; Hparams->projection.sign=+1; break;
			case 2: Hparams->projection.trigo=&Spsi; Hparams->projection.sign=+1; break;
                        //we have the case 3 only with chiral B field
                        case 3: Hparams->projection.trigo=&Noangle; Hparams->projection.sign=+1; break;
	                default: printf("Component is unclear, leaving...\n"); exit(2); break;
		  }
                }
	}
	else
	{
		printf("Nucleus identification is not clear, exiting...\n");
		exit(2);
	}
    
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (intervals_limit);
    gsl_set_error_handler_off();
  
    gsl_function F;
    F.function = &H0_integrand;
    F.params = Hparams;
 
    integration_error=gsl_integration_qags(&F, xl, xr, gsl_abserr, gsl_relerr, intervals_limit, w, &result, &uncertainty); 
    if(integration_error != 0) {
       failed=-1;
       integration_error=gsl_integration_qags(&F, xl, xr, gsl_abserr*coarse_error_factor1, gsl_relerr*coarse_error_factor1, intervals_limit, w, &result, &uncertainty); 
       if(integration_error != 0) {
         failed=-2;
         integration_error=gsl_integration_qags(&F, xl, xr, gsl_abserr*coarse_error_factor2, gsl_relerr*coarse_error_factor2, intervals_limit, w, &result, &uncertainty); 
         if(integration_error != 0) failed=-3;
       }
    }
	
    gsl_integration_workspace_free (w);


    if(Hparams->Btype == 1) {
      return rho*Qbig/(4*M_PI)*v*g*sqrh*result; }
    else {
         if(Hparams->component==3) {
           return schir*rho*Qbig*g*v/(8*M_PI*sqrh)*result;}
         else { //we are computing x or y chiral components
           return (-schir)*rho*Qbig*g*g*v/(8*M_PI*sqrh)*result; }
         }
}
 
double dist1(const double x, const double y)
{ 
    return sqrt(R*R - (x+B/2.)*(x+B/2.)-y*y);
}

double dist2(const double x, const double y)
{ 
	return sqrt(R*R - (x-B/2.)*(x-B/2.)-y*y);
}

double zpos1(const double t, const double z, const double v)
{
	return v*t+z;
}

double zpos2(const double t, const double z, const double v)
{
	return v*t-z;
}


double Noangle(const double x, const double y, const double xi, const double yi)
{
  return 1;
}

double Cpsi(const double x, const double y, const double xi, const double yi)
{
	if((x==xi) && (y==yi))
	{
		return 0;
	}
	else
	{
		return (x-xi)/sqrt((x-xi)*(x-xi)+(y-yi)*(y-yi));
	}
}

double Spsi(const double x, const double y, const double xi, const double yi)
{
	if((x==xi) && (y==yi)) 
	{
		return 0;
	}
	else
	{
		return (y-yi)/sqrt((x-xi)*(x-xi)+(y-yi)*(y-yi));
	}
	
}

double H1_internal_integrand(double yi, void *p)
{
    double result;
	parameters * params = (parameters *)p;
	double t = (params->t);
    double x = (params->x);
    double y = (params->y);
    double z = (params->z);
    double xi= (params->xi);
    params->yi=yi;

    return 2*params->dist(xi,yi)*Hpoint(params)*(params->projection.sign*params->projection.trigo(x,y,xi,yi));
 }   

double H0_integrand(double xi, void *p)
{
    double result;
	parameters * params = (parameters *)p;
	double t = (params->t);
    double x = (params->x);
    double y = (params->y);
    double z = (params->z);
    params->xi=xi;

    return H1_internal(params);
 }   

void print_syntax_and_exit()
{
  printf("Syntax:\n./Hscan.exe b_impact seq Borigin(0,1,2) outputfile z1 [z2] [z3]...\n or \n./Hscan.exe b_impact range Borigin(0,1,2) outputfile zmin zmax Npoints\n or \n./Hscan.exe b_impact point Borigin(0,1,2) x y z\n");
  exit(2);
}
