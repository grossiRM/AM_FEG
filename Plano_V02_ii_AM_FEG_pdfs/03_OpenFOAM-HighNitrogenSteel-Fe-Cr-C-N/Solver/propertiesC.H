		
		Ysat = Y;
	    	
       forAll(mesh.C(), i)
        {
				// DN solid phase (below 100°C no diffusion of nitrogen) 
 
				//if ((T[i] >= 100) && (T[i] < Tl.value()))  
				if ((T[i] >= 500) && (T[i] < Tl.value()))  
					{
						//Info << "Dc ist nicht Null" << endl;
						DNs = Foam::pow(10, -9515/T[i] + 0.4588) * 0.0001 * Dunit; //Nitrogen diffusion in vacancy-rich ferrite and austenite, from first principles to applications
					}
				else
					{
						DNs = VSMALL*Dunit;
					}
				
				if (T[i] > Ts.value()) 
					{
						DNl = (8.866e-12*T[i]-1.214e-08) * Dunit; // Raw Iron Coefficient of diffusion at 1600°C y-Fe Metallkundliche und metallurgische Grundlagen zur Entwicklung stickstoffreicher, zäher, hochfester austenitischer Stähle
					}
				else
					{
						DNl = VSMALL*Dunit;
					}
					
				DN[i] = (DNl * lf[i] + DNs* (1-lf[i])).value();  
				
				if (T[i] > Ts.value() && T[i] <= Tl.value())
				   {
						Ysat[i] = -6.836e-06*T[i] +0.01418;
				   }
				if (T[i] > Tl.value() && T[i] <= 2273.15)
				   {
						Ysat[i] = 0.009778*Foam::exp(-0.0008539*T[i]);
				   }
				if (T[i] > 2273.15)
				   {
						Ysat[i] = 0.0014;
				   }	
				if (T[i] > Ts.value() && T[i] <= T_v.value())
				   {
						k_rate[i] = Sh.value()*DN[i]/Sh_L.value();
				   }	 
				if (T[i] > T_v.value())
				   {
						k_rate[i] = Sh.value()*(8.866e-12*T_v.value()-1.214e-08)/Sh_L.value();
				   }
				if (T[i] <= Ts.value())
				   {
						k_rate[i] = 0;
				   }	    
	    }
	    
	    DN *= pos(alpha1 - 0.01);
	    DN += VSMALL*  Dunit;
	    DN.dimensions().reset(phi.dimensions()/dimLength);
	    lVol = max(lVol, pos0(lf - 0.01));
	    
		Info << "min/max(DN) = " << min(DN).value() << ", " << max(DN).value() << endl;
 		Info << "min/max(Ysat) = " << min(Ysat).value() << ", " << max(Ysat).value() << endl;
 		Info << "min/max(krate) = " << min(k_rate).value() << ", " << max(k_rate).value() << endl;


               
          
          
        
