	#include <math.h>
	#include <stdio.h>
	//#include "engine.h"
	//#include <string.h>
	//#define BUFSIZE 256

	double centralDiffX(double firstX, double secondX, double deltaX){
	return (firstX - secondX)/(2*deltaX);
	}

	double centralDiffY(double firstY, double secondY, double deltaY){
	return (firstY - secondY)/(2*deltaY);
	}

	double forwardbackDiff(double first, double second, double deltaF){
	return (first - second)/(deltaF);
	}

	double secondOrderCentral(double front, double middle, double last, double deltaSecond){
	return (front -(2*middle) + last)/(pow(deltaSecond,2));
	}


	int main(){
	//DECLARE VARIABLES
	int nx = 5; //step in x
	int ny = 5; //step in y
	int nt = 100; //step in time
	int nit = 50; //iterations
	int c = 1; 
	double dx = 2/(double)(nx-1);
	double dy = 2/(double)((ny-1));
	double x[nx]; //x vector
	double y[ny]; //y vector
	double rho = 1; //density
	double nu = 0.1; //viscosity
	double dt = 0.001; //time 
	double u_wall = 1.0;
	double u[nx][ny]; //U component
	double v[nx][ny]; //V component
	double un[nx][ny]; //holder for U
	double vn[nx][ny]; //holder for V
	double b[nx][ny]; //square brackets of Pressure Poisson Equation
	double p[nx][ny]; //Pressure field
	double pn[nx][ny]; //Pressure field


	//Setting zero matrices (all rows, columns)

	for (int x = 0; x<nx;x++){
		for (int y = 0; y<ny;y++){
			u[x][y] = 0;
			v[x][y] = 0;
			b[x][y] = 0;
			p[x][y] = 0;
			pn[x][y] = 0;
			un[x][y] = 0;
			vn[x][y] = 0;
		}
	}

	//Populating X

	for (int i = 0; i<nx; i++){
		x[0] = 0;
		x[i+1] = x[i] + dx;
	}

	//Populating Y

	for (int j = 0; j<ny; j++){
		y[0] = 0;
		y[j+1] = y[j] + dy;
	}

	
	//MAIN LOOP

	for (int it = 0; it<nt; it++){

		//A
		//CALCULATE SQUARE BRACKET SIDE, 'B'
		for (int i = 1; i<nx-1;i++){
			for (int j = 1; j<ny-1;j++){
				b[i][j] = ((rho*-1)*((pow(centralDiffX(u[i+1][j],u[i-1][j],dx),2))+(pow(centralDiffY(v[i][j+1],u[i][j-1],dy),2))+(2*(centralDiffY(u[i][j+1],u[i][j-1],dy))*centralDiffX(v[i+1][j],v[i-1][j],dx)))) + ((rho/dt)*(forwardbackDiff(u[i][j],u[i-1][j],dx))+(forwardbackDiff(v[i][j],v[i][j-1],dy)));
			}
		}	
	


		//B
		//CALCULATE PRESSURE
		for (int iit=0; iit<nit; iit++){
			//1. copy pn to p
			for (int xiter = 0; xiter<nx;xiter++){
				for (int yiter = 0; yiter<ny; yiter++){
					pn[xiter][yiter] = p[xiter][yiter];
				}
			}

			//2. calculate p 

			for (int i = 0; i<nx;i++){
				for (int j = 0; j<ny;j++){
		 			p[i][j] = (b[i][j]-(pow(dy,2)*(pn[i+1][j]+pn[i-1][j]))-(pow(dx,2)*(pn[i][j+1]-pn[i][j-1])))/(-2*((pow(dx,2))+pow(dy,2)));
				}
			}

			//3. boundary condition for p

			for (int yrange = 0; yrange <ny;yrange++){
				p[ny-1][yrange] = p[ny-2][yrange]; //dp/dy = 0 at y = 2
				p[0][yrange] = p[1][yrange]; // dp/dy = 0 at y = 0
				p[ny-1][yrange] = 0;			
			}

			for (int xrange = 0; xrange <nx;xrange++){
				p[xrange][0] = p[xrange][1]; //dp/dx = 0 at x = 0			
			}
		}

		//C

		//CALCULATE U AND V

		//1. Copy un to u, vn to v

		for (int xiteration = 0; xiteration<nx;xiteration++){
			for (int yiteration = 0; yiteration<ny; yiteration++){
				vn[xiteration][yiteration] = v[xiteration][yiteration];
				un[xiteration][yiteration] = u[xiteration][yiteration];
			}
		}

		//2. calculate V field

		for (int i =1; i<nx-1;i++){
			for (int j = 1; j<ny-1; j++){
				u[i][j] = un[i][j]-(un[i][j]*dt*(forwardbackDiff(un[i][j],un[i-1][j],dx)))-(vn[i][j]*dt*(forwardbackDiff(un[i][j],un[i][j-1],dy)))-((dt/rho)*(centralDiffX(p[i+1][j],p[i-1][j],dx)))+((nu*dt)*(secondOrderCentral(un[i+1][j],un[i][j],un[i-1][j],dx))+((secondOrderCentral(un[i][j+1],un[i][j],un[i][j-1],dy))));

				v[i][j] = vn[i][j]-(un[i][j]*dt*(forwardbackDiff(vn[i][j],vn[i-1][j],dx)))-(vn[i][j]*dt*(forwardbackDiff(vn[i][j],vn[i][j-1],dy)))-((dt/rho)*(centralDiffX(p[i][j+1],p[i][j-1],dy)))+((nu*dt)*(secondOrderCentral(vn[i+1][j],vn[i][j],vn[i-1][j],dx))+((secondOrderCentral(vn[i][j+1],vn[i][j],vn[i][j-1],dy))));
			}
		}

		//3. Boundary condition for U, V
		for (int xbound = 0; xbound <nx; xbound++){
			for (int ybound = 0; ybound <ny; ybound++){
				u[0][ybound] = u_wall;
				u[xbound][0] = 0;
				u[nx-1][ybound] = 0;
				u[xbound][ny-2] = 0;


				v[0][ybound] = 0;
				v[ny-1][ybound] = 0;
				v[xbound][0] = 0;
				v[nx-1][ybound] = 0;
			}
		}		
	}	
	 
	
	
	//Testing printing out matrix, still has NAN 

	
		for (int es = 0; es <nx; es++){
			for (int why = 0; why<ny; why++){
				printf("%1.3f\t", u[es][why]);
			}
			printf("\n");
		}
	

	return 0;
	
}
	 



	/*

	//MATLAB PLOTTING

	//Declare basic pointers
	Engine *ep;
	mxArray *p_pointer = NULL;
	
	mxArray *x_pointer = NULL;
	mxArray *u_pointer = NULL;
	mxArray *v_pointer = NULL;

	mxArray *y_pointer = NULL, *result = NULL;
	char buffer[BUFSIZE+1];

	//Test if Matlab is opened or not
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return EXIT_FAILURE;
	}

	x_pointer = mxCreateDoubleMatrix(1,nx,mxREAL);
	memcpy((void *)mxGetPr(x_pointer), (void *)x, sizeof(x));
	engPutVariable(ep, "X_matlab", x_pointer);

	y_pointer = mxCreateDoubleMatrix(1,ny,mxREAL);
	memcpy((void *)mxGetPr(y_pointer), (void *)y, sizeof(y));
	engPutVariable(ep, "Y_matlab", y_pointer);

	p_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(p_pointer), (void *)p, sizeof(p));
	engPutVariable(ep, "P_matlab", p_pointer);

	u_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(u_pointer), (void *)u, sizeof(u));
	engPutVariable(ep, "U_matlab", u_pointer);

	v_pointer = mxCreateDoubleMatrix(nx,ny,mxREAL);
	memcpy((void *)mxGetPr(v_pointer), (void *)v, sizeof(v));
	engPutVariable(ep, "V_matlab", v_pointer);

	//MESH GRID

	engEvalString(ep,"[X Y] = meshgrid(X_matlab,Y_matlab);");
	//engEvalString(ep,"contour(X_matlab,Y_matlab,P_matlab);");
	engEvalString(ep,"quiver(X_matlab,Y_matlab,U_matlab,V_matlab);");
	fgetc(stdin);


	*/  




	















