WITH  Ada.Numerics.Generic_Elementary_Functions,
      Text_IO,
      Provide_partition_function,
      Continuous_Opacity;           -- exports  TYPE  Level_Arr
PROCEDURE  Copcon  IS

  --    Program to calculate wavelength- and depth-dependent continuous
  --    opacities for input in Ada line synthesis programs provided by
  --    M.J. Stift
  --    opacity sources and structure largely taken from
  --    Chmielewski, Publ. Obs. Geneve, Serie B, Fasc. 7 (1979)
  --
  --    Written 1996 by Martin J. Stift
  --                    Institut fuer Astronomie
  --                    Tuerkenschanzstr. 17
  --                    A-1180 Wien
  --                    AUSTRIA
  --
  --                    e-mail:  stift@astro.univie.ac.at
  --
  --    This program is free software; you can redistribute it and/or modify
  --    it under the terms of the GNU General Public License as published by
  --    the Free Software Foundation; either version 2 of the License, or
  --    (at your option) any later version.
  --
  --    This program is distributed in the hope that it will be useful,
  --    but WITHOUT ANY WARRANTY; without even the implied warranty of
  --    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  --    GNU General Public License for more details.
  --
  --    You should have received a copy of the GNU General Public License
  --    along with this program; if not, write to the Free Software
  --    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  No_real_digits : CONSTANT   := 11;

  Ilayer_a      : Integer;
  U_2, U_4, U_5 : Text_IO.File_Type;

  PACKAGE  Int_IO  IS NEW  Text_IO.Integer_IO (Integer);

  -------------------------------------------------------
  PROCEDURE  Copcon_Main (Ilayer_a : IN Integer)  IS

    PRAGMA OPTIMIZE (Time);

    TYPE  Real    IS DIGITS  No_real_digits;

    TYPE  Real_Arr  IS ARRAY (Integer RANGE <>)  OF  Real;

    PACKAGE Real_IO IS NEW Text_IO.FLOAT_IO (Real);
    PACKAGE Math    IS NEW Ada.Numerics.Generic_Elementary_Functions (Real);
    USE     Math;
    PACKAGE Opac    IS NEW Continuous_Opacity(Real);
    PACKAGE Part    IS NEW Provide_partition_function (Real, Real_Arr);

    TYPE  Atom_type  IS RECORD
      Isyb  : Integer;  -- Identification
      Eps_g : Real;  -- Abund., IAU Comm. 12, Grenoble, 1976 (AB(HE) = 0.1).
      Eps_k : Real;  -- Abund.  according to Kurucz (1979)
      Mass  : Real;  -- Atomic masses, Allen, Astrophysical Quantities (1973).
    END RECORD;

    TYPE  Phys_Var_type  IS RECORD
      Tau    : Real;
      Temp   : Real;
      Theta  : Real;
      P_elec : Real;
      P_gas  : Real;
      X      : Real;
      Rho    : Real;
      Rho_x  : Real;
      Kappa  : Real;
    END RECORD;

    TYPE  Concentration_type  IS  RECORD
      HI  : Real;
      HII : Real;
      Hm  : Real;
      H2p : Real;
      H2  : Real;
    END RECORD;

    cLight     : CONSTANT Real  := 2.9979E18;
    cRydberg   : CONSTANT Real  := 13.595;         -- eV
    cBoltzmann : CONSTANT Real  := 1.380622E-16;   -- erg K^-1
    cAtMass    : CONSTANT Real  := 1.66019E-24;
    cStefan    : CONSTANT Real  := 7.564E-15;
    Thomson    : CONSTANT Real  := 66.52;

    No_elements : CONSTANT Integer := Part.No_elements;
    Nma_max     : CONSTANT Integer := Part.Nma_max;
    Iwave       : CONSTANT Integer := 55;
    Idepth      : CONSTANT Integer := 501;
    Ilayer      : CONSTANT Integer := Ilayer_a;

    TYPE  Press_table_type  IS ARRAY (1..Idepth)  OF  Real;

    Lambda : CONSTANT ARRAY (1..Iwave)  OF  Real :=
      (3700.0, 3750.0, 3800.0, 3850.0, 3900.0, 3950.0, 4000.0, 4050.0,
       4100.0, 4150.0, 4200.0, 4250.0, 4300.0, 4350.0, 4400.0, 4450.0,
       4500.0, 4550.0, 4600.0, 4650.0, 4700.0, 4750.0, 4800.0, 4850.0,
       4900.0, 4950.0, 5000.0, 5050.0, 5100.0, 5150.0, 5200.0, 5300.0,
       5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 6000.0, 6100.0,
       6200.0, 6300.0, 6400.0, 6500.0, 6600.0, 6700.0, 6800.0, 6900.0,
       7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0,10000.0         );

    Atmos        : ARRAY (1..Ilayer)  OF  Phys_Var_type;
    Result, Phys : ARRAY (1..Idepth)  OF  Phys_Var_type;

    Abun, Zeta  : Real_Arr(1..No_elements);
    QQ, Chi_ion : Real_Arr(1..Nma_max);

    P_e_table, P_g_table  : Press_table_type;

    X_new, Tau_new, Kappa_new, Temp_new, Rho_new,
    Rhox_new, Pgas_new, Pelec_new, Tab_Theta,
    Tab_Theta_1, Tab_Pelec, Tab_Pelec_1, Tab_Pgas,
    Tab_Pgas_1, Tab_Kappa, Tab_Kappa_1, Tab_Rho,
    Tab_Rho_1, Tab_Rhox, Tab_Rhox_1, Tab_X, Tab_X_1  : Press_table_type := (OTHERS => 0.0);

    Conc      : Concentration_type;
    Departure : Opac.Level_Arr      := (OTHERS => 1.0);

    Flag             : BOOLEAN   := TRUE;
    Stiffness, Ecart : Real      := 1.0;
    Sab, Amu         : Real      := 0.0;

    Rhox_min, Rhox_max, Kap, Density, Freq,
    P_e, V_turb, Alkl, Salph, Cro, Log_g,
    He, T_eff, Help, H_Tau, H_Rho, H_Kappa,
    H_Temp, H_Pgas, H_Pelec                 : Real;
    J_C, J_Mg, J_Al, J_Si, J_Fe,
    In_Last, Ou_Last                        : Integer;
    Pflag, Mat, M_dep                       : Integer         := 1;
    Input_File_name, Output_File_name       : String(1..80);
    Option                                  : CHARACTER;

    --------------------------------------------
    FUNCTION  Log10 (X : Real)  RETURN  Real  IS
    BEGIN
      RETURN  Log (X, 10.0);
    END  Log10;

    -------------------------------------------------------------
    FUNCTION  Indic (X0 : Real;
                     X  :  Press_table_type;
                     M  :  Integer         )  RETURN  Integer  IS

      --  Search for the position of a value X0 in a table of M monotonically
      --  increasing or decreasing values of X.

      --  If K = Indic (X0,X,M)     we have X0 between X(K) and X(K+1) or X0 = X(K).
      --  If X0 is not in the table we have Indic = 1   if X0 is on the side of X(1)
      --                                    Indic = M-1 if X0 is on the side of X(M)

      PRAGMA OPTIMIZE (Time);

      I  : Integer;
      K  : Integer   := 1;
      N  : Integer   := M;
      Dx : Real      := X(N) - X(1);
    BEGIN
      LOOP
        IF  ((N-K-1) = 0)  THEN
          RETURN  K;
        ELSE
          I := (K + N) / 2;
          IF  ((X(I)-X0)*Dx > 0.0)  THEN
            N := I;
          ELSE
            K := I;
          END IF;
        END IF;
      END LOOP;
    END  Indic;

    ----------------------------------------------------------
    FUNCTION  Ypol (X0 : Real;
                    X  : Press_table_type;
                    Y  : Press_table_type;
                    N  : Integer           )  RETURN  Real  IS

      PRAGMA OPTIMIZE (Time);

      S, T, U, V, W : Real;
      K             : Integer;

      --  Quadratic interpolation for the value X0 in the table (X(I),Y(I),I=1,N)
      --  which is monotonically increasing or decreasing in X.
      --  (the extrapolation is linear).

    BEGIN
      K := Indic (X0, X, N);

      --  Linear interpolation
      IF  ((N = 2)  OR  ((X(N)-X0)*(X(1)-X0) >= 0.0))  THEN
        RETURN  Y(K+1) + (Y(K+1)-Y(K)) * (X0-X(K+1)) / (X(K+1)-X(K));

      --  Quadratic interpolation
      ELSE
        K := Integer'Min (K, (N-2));
        U := X0 - X(K);
        V := X0 - X(K+1);
        W := X0 - X(K+2);
        S := (Y(K)   - Y(K+1)) / (X(K)   - X(K+1));
        T := (Y(K+1) - Y(K+2)) / (X(K+1) - X(K+2));
        IF  ((S*T) <= 0.0)  THEN
          RETURN  Y(K) + U * (S + V * (S-T) / (W-U));
        ELSE
          RETURN  Y(K) + U * S / (1.0 - V * (S-T) / (Y(K) - Y(K+2)));
        END IF;
      END IF;
    END  Ypol;

    -----------------------------------------------
    FUNCTION  Sign (X, Y : Real )  RETURN  Real  IS
    BEGIN
      IF  (Y < 0.0)  THEN
        RETURN  -X;
      ELSE
        RETURN   X;
      END IF;
    END;

    ---------------------------------------------
    FUNCTION  Exp10 (X : Real )  RETURN  Real  IS
      Ln10 : CONSTANT Real  := Log (10.0);
    BEGIN
      RETURN  Exp (x * Ln10);
    END;

    ------------------------------
    PROCEDURE  Equation_of_State
     (Temperature : IN   Real;
      Elec_press  : IN   Real;
      Gas_press   : OUT  Real;
      Rho         : OUT  Real)  IS

      --  Physical and thermodynamical variables of a stellar atmosphere.
      --  We are considering the ionisation equilibrium of the metals.
      --  The partition function varies with temperature and electron pressure.
      --  For hydrogen we consider ionisation and the formation of H-, H2, ET H2+ .
      --  See Mihalas, Methods in computational Physics, Vol. 7, p. 8 (1967).

      --  Input data:  Temperature and electronic pressure
      --               Chemical composition. AB(J) gives the abundances of the NEL
      --               elements in Atoms(J).Isyb. cAtMass ist the atomic mass unit.

      --  Results : Gas pressure and density at each depth point
      --            Concentrations, per cm3, of H I, H II, H- (Hm), H2 and H2+ (H2p).
      --            At each depth point I : Zeta(I,J) = N0/U0 = Number of neutral
      --            atoms per cm3 of the element J, divided by the partition function
      --            of the netural metal.
      --            Amen(J) = Abundance of the neutral metal J per neutral
      --            hydrogen for a given point

      --   Ionisation equilibrium of the metals

      PRAGMA OPTIMIZE (Time);

      Theta : Real  := 5040.4 / Temperature;
      DChi  : Real  := 4.98E-04 * Theta * Sqrt (Elec_press);
      Saha  : Real  := 9.0799 - Log10 (Elec_press) - 1.0857 * Log (Theta);
      Sanf  : Real  := 0.0;

      Nma                        : Integer;
      A, B, C, D, E, Pn, Sn,
      F_e, Sfi, HT, FJN, SNFJN,
      c_H, c_H2, c_H2p, c_Hm,
      G2, G3, F1, F2, F3, F4, F5 : Real;

      Amen  : ARRAY (1..No_elements)  OF  Real;

      -------------------------------------------------
      FUNCTION  Root (A, B, C : Real)  RETURN  Real  IS

        PRAGMA OPTIMIZE (Time);

        --  Physical root of a second order equation

        Eps     : Real    := 1.0E-12;
        U, V, D : Real;
      BEGIN
        U := 4.0 * A * C;
        V := B**2;
        IF  (ABS (U / V) - Eps > 0.0)  THEN
          D := Sign (1.0, A) * Sqrt (V - U) - B;
          RETURN  D / (2.0 * A);
        ELSE
          RETURN  - C / B;
        END IF;
      END  Root;

    BEGIN
      FOR  J  IN  2..No_elements  LOOP
        Part.Partition_Function (Iflag   => Pflag,
                                 Iel     => Part.Atoms(J).Isyb,
                                 DChi    => Dchi,
                                 Theta   => Theta,
                                 Q       => QQ,
                                 Chi_ion => Chi_ion,
                                 Nma     => Nma               );
        Pn    := 0.0;
        Sn    := 0.0;
        SNFJN := 0.0;
        FOR  N  IN  1..Nma  LOOP
          FJN   := QQ(N) * Exp10 (Pn);
          SNFJN := SNFJN + (Real (N) - 1.0) * FJN;
          Sn    := Sn + FJN;
          Pn    := Pn + Saha + (Real (N) * DChi - Chi_ion(N)) * Theta;
        END LOOP;
        Zeta(J) := 1.0 / Sn;
        Amen(J) := QQ(1) / Sn;
        Sanf    := Sanf + Abun(J) * SNFJN / Sn;
      END LOOP;

      --  Ionisation and dissociation equilibrium of hydrogen

      Part.Partition_Function (Iflag   => Pflag,
                               Iel     => Part.Atoms(1).Isyb,
                               DChi    => Dchi,
                               Theta   => Theta,
                               Q       => QQ,
                               Chi_ion => Chi_ion,
                               Nma     => Nma         );

      --  equations (27) - (29), (30) for c_H is modified

      c_H   := Exp10 (Saha + (DChi - Chi_ion(1)) * Theta) *
                     QQ(2) / QQ(1) * Elec_press;
      c_H2  := Exp10 (12.533505 - Theta * (4.9251644 - Theta *
                     (0.056191273 - Theta * 0.0032687661)));
      c_H2p := Exp10 (11.206998 - Theta * (2.7942767 + Theta *
                     (0.079196803 - Theta * 0.024790744)));
      c_Hm  := Exp10 (-0.747 * Theta + 2.5 * Log10 (Temperature) + 0.1249);

      G2 := c_H / Elec_press;
      G3 := Elec_press / c_Hm;
      E  := G2 * c_H2 / c_H2p;
      D  := G2 - G3;
      C  := Elec_press / c_H2;
      B  := 2.0 * (1.0 + E);
      A  := 1.0 + G2 + G3;
      F1 := Root (A => C * B**2 + A * D * B - E * A**2,
                  B => 2.0 * A * E - D * B + A * B * Sanf,
                  C => -E - B * Sanf);
      F2 := G2 * F1;
      F3 := G3 * F1;
      F5 := (1.0 - A * F1) / B;
      F4 := E * F5;

      --  equations (35) and (36)
      --  HT denotes number density of hydrogen if it were all in atomic form

      F_e       := F2 - F3 + F4 + Sanf;
      Sfi       := F1 + F2 + F3 + F4 + F5;
      Gas_press := Exp (Log (Elec_press) + Log (1.0 + (Sfi + Sab) / F_e));
      HT        := Elec_press / F_e / (cBoltzmann * Temperature);
      Rho       := cAtMass * Amu * HT;

      Conc.HI  := F1 * HT;
      Conc.HII := F2 * HT;
      Conc.Hm  := F3 * HT;
      Conc.H2p := F4 * HT;
      Conc.H2  := F5 * HT;

      FOR  J  IN  2..No_elements  LOOP
        Zeta(J) := Zeta(J) * Abun(J) * HT;
        Amen(J) := Amen(J) * Abun(J) / F1;
      END LOOP;

    END  Equation_of_State;

    -------------------------------------------------
    PROCEDURE  Model (Input_File_name  : String;
                      Output_File_name : String )  IS

      Dum_1, Dum_2 : Real;
      Read_string  : String (1..4);
      I_Teff       : Integer;
    BEGIN

      --  Read title, number of depth points, T_eff, Log g, helium abundance
      --  and turbulence

      Text_IO.Open  (U_5, Text_IO.In_File,  Input_File_name);
      Text_IO.Create(U_4, Text_IO.Out_File, Output_File_name);

      Text_IO.Set_Col  (U_5, 5);
      Int_IO.Get       (U_5, I_Teff, WIDTH => 7);
      Text_IO.Set_Col  (U_5, 22);
      Real_IO.Get      (U_5, Log_g);
      Text_IO.Skip_Line(U_5);

      T_eff := Real (I_Teff);
      He    := 0.1;

      LOOP
        Text_IO.Get (U_5, Read_string);
        EXIT WHEN  (Read_string = "READ");
        Text_IO.Skip_Line(U_5);
      END LOOP;
      Text_IO.Set_Col  (U_5, 11);
      Int_IO.Get       (U_5, M_dep);
      Text_IO.Skip_Line(U_5);

      Salph := Sab + Abun(1);
      Cro   := cAtMass * Amu / Salph;
      FOR  I  IN  1..M_dep  LOOP
        Real_IO.Get (U_5, Atmos(I).Rho_x);
        Real_IO.Get (U_5, Atmos(I).Temp);
        Real_IO.Get (U_5, Atmos(I).P_gas);
        Real_IO.Get (U_5, Atmos(I).P_elec);
        Real_IO.Get (U_5, Dum_1);
        Real_IO.Get (U_5, Dum_2);
        Real_IO.Get (U_5, V_turb);
        Atmos(I).P_elec := Atmos(I).P_elec * Atmos(I).Temp * cBoltzmann;
        Text_IO.Skip_Line (U_5);
      END LOOP;
      V_turb := V_turb * 1.0e-5;
    END  Model;

    ----------------------------------------------------
    FUNCTION  Opacity (Freq : Real;
                       Temp : Real;
                       P_e  : Real )  RETURN  Real  IS

      --  Continuous opacities at the frequency Freq, at a given point of the
      --  atmosphere with temperature Temp and electron pressure P_e.
      --  Ecart is the departure coefficient from LTE of the ground level of hydrogen
      --  Opacities due to the different forms of hydrogen are calculated according
      --  to the routines of Carbon and Gingerich in Gingerich, ed. (1969) "Theory
      --  and observation of normal stellar atmospheres" (Cambridge, MIT Press).

      --  For the photoionisation of hydrogen the departure from LTE of the
      --  ground level population is taken into account following Gingerich et al.,
      --  Solar Physics 18, 347 (1971).
      --  Opacities due to photoionisation of neutral metals are given by
      --  Vernazza, Avrett and Loeser, ApJS, 30, 1 (1976).
      --  Result : Kappa gives the opacity per cm3 (in cm-1)

      PRAGMA OPTIMIZE (Time);

      USE  Opac;
      Pphi   : Real   := Conc.HI * cBoltzmann * Temp;
      ABSCO  : Real   := 0.0;
      Wave   : Real   := cLight / Freq;
      Theta  : Real   := 5040.4 / Temp;
    BEGIN
      ABSCO := ABSCO + Opac_Hydrogen       (Theta, Freq, Ecart);
      ABSCO := ABSCO + Opac_H_minus        (Theta, Freq) * P_e;
      ABSCO := ABSCO + Opac_Lyman          (Theta, Freq);
      ABSCO := ABSCO + Opac_H2p            (Theta, Freq) * Conc.HII;
      ABSCO := ABSCO + Opac_Quasi_Hydrogen (Theta, Freq) * Conc.HI;
      ABSCO := ABSCO + Opac_Rayleigh_H            (Freq);
      ABSCO := ABSCO + Opac_Rayleigh_H2           (Freq) * Conc.H2 / Conc.HI;
      ABSCO := ABSCO + Chi_C  (Temp, Wave, Zeta(J_C ), Departure)  / Conc.HI;
      ABSCO := ABSCO + Chi_Mg (Temp, Wave, Zeta(J_Mg), Departure)  / Conc.HI;
      ABSCO := ABSCO + Chi_Al (Temp, Wave, Zeta(J_Al), Departure)  / Conc.HI;
      ABSCO := ABSCO + Chi_Si (Temp, Wave, Zeta(J_Si), Departure)  / Conc.HI;
      ABSCO := ABSCO + Chi_Fe (Temp, Wave ,Zeta(J_Fe), Departure)  / Conc.HI;
      ABSCO := ABSCO + Thomson * P_e / Pphi;

      --  KAPPA per gram = OPV / RHO
      --  KAPPA per atom = OPV * cBoltzmann * Temp / (P_g - P_e)

      RETURN  ABSCO * Conc.HI * 1.0E-26;
    END  Opacity;

    --------------------------------------------------------
    PROCEDURE  Electron_pressure (Temp    : IN   Real;
                                  Pr_gas  : IN   Real;
                                  Pr_ele  : OUT  Real;
                                  Kappa   : OUT  Real;
                                  Rho     : OUT  Real  )  IS

      --  Determine electronic pressure corresponding to a given gas pressure
      --  Then calculcate standard opacity (per gram at 5000 A).

      PRAGMA OPTIMIZE (Time);

      Kappa_5000, Density, P_g : Real;
      Freq                     : Real    := cLight / 5000.0;
      M                        : Integer := 0;

      ----------------------------------
      PROCEDURE  Press_Order
       (M       : IN OUT  Integer;
        G_press : IN      Real;
        E_press : IN      Real     )  IS

        PRAGMA OPTIMIZE (Time);

        iflag : BOOLEAN;
      BEGIN
        iflag := TRUE;
        IF  M = 0 THEN
          iflag      := FALSE;
          P_e_table(M+1) := E_press;
          P_g_table(M+1) := G_press;
          M              := M + 1;
        ELSE
          FOR  I1  IN  1..M  LOOP
            IF  (G_press < P_g_table(I1))  THEN
              iflag := FALSE;
              FOR  I2  IN REVERSE  I1..M  LOOP
                P_e_table(I2+1) := P_e_table(I2);
                P_g_table(I2+1) := P_g_table(I2);
              END LOOP;
              P_e_table(I1) := E_press;
              P_g_table(I1) := G_press;
              M             := M + 1;
              EXIT;
            END IF;
          END LOOP;
        END IF;
        IF  iflag  THEN
          M            := M + 1;
          P_e_table(M) := E_press;
          P_g_table(M) := G_press;
        END IF;
      END  Press_Order;

    BEGIN
      P_e := Pr_gas * 1.0E-3;
      FOR  I  IN  1..3  LOOP
        Equation_of_State (Temperature => Temp,
                           Elec_press  => P_e,
                           Gas_press   => P_g,
                           Rho         => Density  );
        Press_Order (M       => M,
                     G_press => Log (P_g),
                     E_press => Log (P_e)    );
        P_e := 2.0 * P_e;
      END LOOP;
      LOOP
        P_e :=  Exp (Ypol (X0 => Log (Pr_gas),
                           X  => P_g_table,
                           Y  => P_e_table,
                           N  => M           ) );
        Equation_of_State (Temperature => Temp,
                           Elec_press  => P_e,
                           Gas_press   => P_g,
                           Rho         => Density);
        Press_Order (M       => M,
                     G_press => Log (P_g),
                     E_press => Log (P_e));
        EXIT WHEN  (ABS ((P_g - Pr_gas) / Pr_gas) < 1.0E-4);
      END LOOP;
      Kappa_5000 := Opacity (Freq, Temp, P_e) / Density;
      Kappa      := Kappa_5000;
      Rho        := Density;
      Pr_ele     := P_e;
    END  Electron_pressure;

  BEGIN

    Text_IO.Skip_Line(U_2); Text_IO.Set_Col (U_2, 35);
    Text_IO.Get_Line (U_2, ITEM => Input_File_name,  LAST => In_Last);
    Text_IO.Set_Col  (U_2, 35);
    Text_IO.Get_Line (U_2, ITEM => Output_File_name, LAST => Ou_Last);
    Text_IO.Set_Col  (U_2, 35);
    Text_IO.Get      (U_2, Option);
    Text_IO.Skip_Line(U_2); Text_IO.Set_Col (U_2, 35);
    Int_IO.Get       (U_2, Mat);
    Text_IO.Skip_Line(U_2); Text_IO.Set_Col (U_2, 35);
    Real_IO.Get      (U_2, Rhox_min);
    Text_IO.Skip_Line(U_2); Text_IO.Set_Col (U_2, 35);
    Real_IO.Get      (U_2, Rhox_max);

    -- interpolate Mat points in rho scale

    Rhox_min := Log10 (Rhox_min);
    Rhox_max := Log10 (Rhox_max);
    FOR  K  IN  1..Mat  LOOP
      Rhox_new(K) := Rhox_min + Real (K-1) * (Rhox_max - Rhox_min) / Real (Mat-1);
    END LOOP;

    -- default : IAU abundances

    IF  ((Option = 'K')  OR  (Option = 'k'))  THEN
      FOR  J  IN  1..No_elements  LOOP
        Abun(J) := Exp10 (Part.Atoms(J).Eps_k - 12.0);
      END LOOP;
    ELSE
      FOR  J  IN  1..No_elements  LOOP
        Abun(J) := Exp10 (Part.Atoms(J).Eps_g - 12.0);
      END LOOP;
    END IF;

    Amu := Abun(1) * Part.Atoms(1).Mass;
    FOR  J  IN  2..No_elements  LOOP
      Sab := Sab + Abun(J);
      Amu := Amu + Abun(J) * Part.Atoms(J).Mass;
      IF  (Part.Atoms(J).Isyb =  6)  THEN  J_C  := J;  END IF;
      IF  (Part.Atoms(J).Isyb = 12)  THEN  J_Mg := J;  END IF;
      IF  (Part.Atoms(J).Isyb = 13)  THEN  J_Al := J;  END IF;
      IF  (Part.Atoms(J).Isyb = 14)  THEN  J_Si := J;  END IF;
      IF  (Part.Atoms(J).Isyb = 26)  THEN  J_Fe := J;  END IF;
    END LOOP;
    Amu := Amu / Part.Atoms(1).Mass;

    -- read input model

    Model (Input_File_name(1..In_Last),
           Output_File_name(1..Ou_Last) );

    -- determine electron pressure for given gas pressure

    FOR  K  IN  1..M_dep  LOOP
      Electron_pressure (Temp   => Atmos(K).Temp,
                         Pr_gas => Atmos(K).P_gas,
                         Pr_ele => Result(K).P_elec,
                         Kappa  => Result(K).Kappa,
                         Rho    => Result(K).Rho    );
      Int_IO.Put  (K);
      Text_IO.Put ("    ");
      Real_IO.Put ((Atmos(K).P_elec-Result(K).P_elec)/Atmos(K).P_elec);
      Text_IO.Put ("    ");
      Text_IO.New_Line;
    END LOOP;

    Int_IO.Put  (U_4, Mat, WIDTH => 4);
    Text_IO.New_Line (U_4);
    Real_IO.Put (U_4, T_eff,  FORE => 8, AFT => 1, EXP => 0);
    Real_IO.Put (U_4, Log_g,  FORE => 8, AFT => 2, EXP => 0);
    Real_IO.Put (U_4, He,     FORE => 4, AFT => 2, EXP => 0);
    Text_IO.New_Line (U_4);

    FOR  K  IN  1..M_dep  LOOP
      Tab_Kappa(K) := Result(K).Kappa;
      Tab_Theta(K) := Atmos(K).Temp;
      Tab_Pelec(K) := Atmos(K).P_elec;
      Tab_Pgas(K)  := Atmos(K).P_gas;
      Tab_Rho(K)   := Result(K).Rho;
      Tab_Rhox(K)  := Atmos(K).Rho_x;
    END LOOP;
    Text_IO.Put_Line(U_4, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

    -- interpolate

    FOR  K  IN  1..Idepth  LOOP
      Rhox_new(K) := Exp10 (Rhox_min + Real (K-1) *  (Rhox_max - Rhox_min) /
                     Real (Idepth-1));
    END LOOP;
    FOR  K  IN  1..Idepth  LOOP
      Rho_new(K)   := Ypol (Rhox_new(K), Tab_Rhox, Tab_Rho,   M_dep);
      Kappa_new(K) := Ypol (Rhox_new(K), Tab_Rhox, Tab_Kappa, M_dep);
      Temp_new(K)  := Ypol (Rhox_new(K), Tab_Rhox, Tab_Theta, M_dep);
      Pgas_new(K)  := Ypol (Rhox_new(K), Tab_Rhox, Tab_Pgas,  M_dep);
      Pelec_new(K) := Ypol (Rhox_new(K), Tab_Rhox, Tab_Pelec, M_dep);
    END LOOP;
    Tau_new(1) := 0.0;
    FOR  K  IN  2..Idepth  LOOP
      Tau_new(K) := (Rhox_new(K)  - Rhox_new(K-1)) *
                    (Kappa_new(K) + Kappa_new(K-1)) / 2.0;
    END LOOP;

    FOR  K  IN  2..Idepth LOOP
      Tau_new(K) := Tau_new(K) + Tau_new(K-1);
    END LOOP;

    -- output interpolated model atmosphere with kappa_5000

    Rhox_min := Rhox_min + (Rhox_max - Rhox_min) / 500.0;
    FOR  K  IN  1..Mat  LOOP
      Help    := Exp10 (Rhox_min + Real (K-1) *  (Rhox_max - Rhox_min) /
                          Real (Mat-1));
      H_Tau   := Ypol (Help, Rhox_new, Tau_new,   Idepth);
      H_Rho   := Ypol (Help, Rhox_new, Rho_new,   Idepth);
      H_Kappa := Ypol (Help, Rhox_new, Kappa_new, Idepth);
      H_Temp  := Ypol (Help, Rhox_new, Temp_new,  Idepth);
      H_Pgas  := Ypol (Help, Rhox_new, Pgas_new,  Idepth);
      H_Pelec := Ypol (Help, Rhox_new, Pelec_new, Idepth);
      Alkl    := -Log10 ((H_Temp * H_Kappa * H_Rho) * cBoltzmann /
                         (H_Pgas - H_Pelec));
      Real_IO.Put (U_4, H_Tau,           FORE => 4, AFT => 3, EXP => 3);
      Real_IO.Put (U_4, 5040.4/H_Temp,   FORE => 4, AFT => 4, EXP => 0);
      Real_IO.Put (U_4, Log10 (H_Pelec), FORE => 4, AFT => 3, EXP => 0);
      Real_IO.Put (U_4, Log10 (H_Pgas),  FORE => 4, AFT => 3, EXP => 0);
      Real_IO.Put (U_4, Alkl,            FORE => 4, AFT => 3, EXP => 0);
      Text_IO.New_Line (U_4);
      Phys(K).Temp   := H_Temp;
      Phys(K).P_gas  := H_Pgas;
      Phys(K).P_elec := H_Pelec;
    END LOOP;

    -- output relative error between input P_elec and calculated P_elec

    FOR  K  IN  1..Mat  LOOP
      Electron_pressure (Temp   => Phys(K).Temp,
                         Pr_gas => Phys(K).P_gas,
                         Pr_ele => H_Pelec,
                         Kappa  => Phys(K).Kappa,
                         Rho    => Phys(K).Rho   );
      Phys(K).Theta := 5040.4 / Phys(K).Temp;
      Int_IO.Put (K);
      Text_IO.Put ("    ");
      Real_IO.Put ((Phys(K).P_elec - H_Pelec) / Phys(K).P_elec);
      Text_IO.Put ("    ");
      Text_IO.New_Line;
    END LOOP;

    -- for all wavelengths, calculate and output continuous opacities

    Int_IO.Put (U_4, Iwave, WIDTH => 3);
    Text_IO.New_Line (U_4);
    FOR  L  IN  Lambda'RANGE  LOOP
      Freq := cLight / Lambda(L);
      Real_IO.Put (U_4, Lambda(L),  FORE => 7, AFT => 1, EXP => 0);
      Text_IO.New_Line (U_4);
      FOR  K  IN  1..Mat  LOOP
        Equation_of_State (Temperature => Phys(K).Temp,
                           Elec_press  => Phys(K).P_elec,
                           Gas_press   => Help,
                           Rho         => Density  );
        Kap   := Opacity (Freq, Phys(K).Temp, Phys(K).P_elec);
        Alkl  := -Log10 (Phys(K).Temp * Kap * cBoltzmann /
                         (Phys(K).P_gas - Phys(K).P_elec));
        Real_IO.Put (U_4, Alkl, FORE => 3, AFT => 3, EXP => 0);
        IF  (((K MOD 10) = 0)  OR  (K = Mat))  THEN
          Text_IO.New_Line (U_4);
        END IF;
      END LOOP;
    END LOOP;

  END  Copcon_Main;

BEGIN
   Text_IO.Open   (U_2, Text_IO.In_File,  "cop_in.dat");

   Text_IO.Skip_Line (U_2);
   Text_IO.Skip_Line (U_2);
   Text_IO.Skip_Line (U_2); Text_IO.Set_Col (U_2, 35);
   Int_IO.Get        (U_2, Ilayer_a);  -- max. no. original depth points
                                       -- from model atmosphere
                                       -- actual number will be M_dep
   Copcon_Main (Ilayer_a);

END  Copcon;
GENERIC
  TYPE  Real  IS DIGITS <>;
PACKAGE  Continuous_Opacity  IS

  --    Continuous opacity package (specification part)
  --
  --    Written 1996 by Martin J. Stift
  --                    Institut fuer Astronomie
  --                    Tuerkenschanzstr. 17
  --                    A-1180 Wien
  --                    AUSTRIA
  --
  --                    e-mail:  stift@astro.univie.ac.at
  --
  --    This program is free software; you can redistribute it and/or modify
  --    it under the terms of the GNU General Public License as published by
  --    the Free Software Foundation; either version 2 of the License, or
  --    (at your option) any later version.
  --
  --    This program is distributed in the hope that it will be useful,
  --    but WITHOUT ANY WARRANTY; without even the implied warranty of
  --    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  --    GNU General Public License for more details.
  --
  --    You should have received a copy of the GNU General Public License
  --    along with this program; if not, write to the Free Software
  --    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  TYPE  Const_Arr  IS PRIVATE;
  TYPE  Level_Arr  IS ARRAY (1..8)   OF  Real;

  -------------------------------------------------
  FUNCTION  Para (X     : Const_Arr;
                  Y     : Const_Arr;
                  Arg_x : Real;
                  N     : Integer  )  RETURN  Real;

  ----------------------------------------------------
  FUNCTION  Gaunt (Frequency : Real;
                   Level     : Integer)  RETURN  Real;

  ---------------------------------------------------------
  FUNCTION  Opac_Hydrogen (Theta     : Real;
                           Frequency : Real;
                           Ecart     : Real)  RETURN  Real;

  ----------------------------------------------------
  FUNCTION  Opac_H2p (Theta     : Real;
                      Frequency : Real)  RETURN  Real;

  --------------------------------------------------------
  FUNCTION  Opac_H_minus (Theta     : Real;
                          Frequency : Real)  RETURN  Real;

  ------------------------------------------------------
  FUNCTION  Opac_Lyman (Theta     : Real;
                        Frequency : Real)  RETURN  Real;

  -----------------------------------------------------------
  FUNCTION  Opac_Rayleigh_H (Frequency : Real)  RETURN  Real;

  ------------------------------------------------------------
  FUNCTION  Opac_Rayleigh_H2 (Frequency : Real)  RETURN  Real;

  ---------------------------------------------------------------
  FUNCTION  Opac_Quasi_Hydrogen (Theta     : Real;
                                 Frequency : Real)  RETURN  Real;

  ---------------------------------------------------
  FUNCTION  Chi_C (Temp   : Real;
                   Wave   : Real;
                   Zeta_0 : Real;
                   Bl     : Level_Arr)  RETURN  Real;

  ----------------------------------------------------
  FUNCTION  Chi_Mg (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real;

  ----------------------------------------------------
  FUNCTION  Chi_Al (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real;

  ----------------------------------------------------
  FUNCTION  Chi_Si (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real;

  ----------------------------------------------------
  FUNCTION  Chi_Fe (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real;

PRIVATE
  TYPE  Const_Arr  IS ARRAY (1..46)  OF  Real;
END  Continuous_Opacity;
WITH  Ada.Numerics.Generic_Elementary_Functions;
PACKAGE BODY  Continuous_Opacity  IS

  --    Continuous opacity package (body)
  --
  --    Written 1996 by Martin J. Stift
  --                    Institut fuer Astronomie
  --                    Tuerkenschanzstr. 17
  --                    A-1180 Wien
  --                    AUSTRIA
  --
  --                    e-mail:  stift@astro.univie.ac.at
  --
  --    This program is free software; you can redistribute it and/or modify
  --    it under the terms of the GNU General Public License as published by
  --    the Free Software Foundation; either version 2 of the License, or
  --    (at your option) any later version.
  --
  --    This program is distributed in the hope that it will be useful,
  --    but WITHOUT ANY WARRANTY; without even the implied warranty of
  --    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  --    GNU General Public License for more details.
  --
  --    You should have received a copy of the GNU General Public License
  --    along with this program; if not, write to the Free Software
  --    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  PACKAGE  Math  IS NEW  Ada.Numerics.Generic_Elementary_Functions (Real);
  USE      Math;

  cLight    : CONSTANT  Real  := 2.9979E18;
  C2        : CONSTANT  Real  := 1.43883E8;
  Ecart     : Real            := 1.0;
  Departure : Level_Arr       := (OTHERS => 1.0);

  ------------------------------------------------------
  FUNCTION  Para (X     : Const_Arr;
                  Y     : Const_Arr;
                  Arg_x : Real;
                  N     : Integer    )  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Parabolic interpolation - Gingerich
    --  Gives para, the interpolated value in the X array for argument Arg_x.
    --  For 0.0 Jx, the A, B, C constants will be unchanged.

    A, B, C       : Real;
    X12, X23, X31 : Real;
    Y12           : Real;
  BEGIN
    X12 := X(N-1) - X(N);
    X23 := X(N)   - X(N+1);
    X31 := X(N+1) - X(N-1);
    Y12 := Y(N-1) - Y(N);
    C   := (X12 * (Y(N) - Y(N+1)) - X23 * Y12) / (X12 * X23 * X31);
    B   := Y12 / X12 - C * (X(N-1) + X(N));
    A   := Y(N-1) - (B + C * X(N-1)) * X(N-1);
    RETURN  A + (B + C * Arg_x) * Arg_x;
  END  Para;

  -------------------------------------------------------
  FUNCTION  Gaunt (Frequency : Real;
                   Level     : Integer)  RETURN  Real  IS

    --  Gaunt factors for b-f hydrogen, Gingerich, March 1964

    PRAGMA OPTIMIZE (Time);

    W : Real  := cLight * 1.0E-3 / Frequency;
  BEGIN
    CASE  Level  IS
      WHEN  1  =>
        RETURN  0.9916 + ( 9.068E-3 - 0.2524   * W) * W;
      WHEN  2  =>
        RETURN  1.105  + (-7.922E-2 + 4.536E-3 * W) * W;
      WHEN  3  =>
        RETURN  1.101  + (-3.290E-2 + 1.152E-3 * W) * W;
      WHEN  4  =>
        RETURN  1.101  + (-1.923E-2 + 5.110E-4 * W) * W;
      WHEN  5  =>
        RETURN  1.102  + (-0.01304  + 2.638E-4 * W) * W;
      WHEN  6  =>
        RETURN  1.0986 + (-0.00902  + 1.367E-4 * W) * W;
      WHEN  OTHERS =>
        RETURN  1.0;
    END CASE;
  END  Gaunt;

  ------------------------------------------------------------
  FUNCTION  Opac_Hydrogen (Theta     : Real;
                           Frequency : Real;
                           Ecart     : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Hydrogen opacity per neutral atom, scaled by E26.
    --  Based on Menzel and Pekeris, MNRAS 96, 77, 1935.
    --  Gingerich - 21 July 1960 - 21 March 1967.

    T       : Real   := 5040.4 / Theta;
    U       : Real   := 157779.0 / T;
    Ehvkt   : Real   := Exp (-4.79895E-11 / T * Frequency);
    D       : Real   := 0.0;
    W       : Real   := cLight / Frequency;
    Square  : Real;
    B, G3   : Real;
    N       : Integer  := 0;
    N_prime : Integer;
  BEGIN
    IF  (Frequency < 0.5137538E14)  THEN
      N := 9;

    --  determine lowest atomic level N considered
    ELSE
      LOOP
        N := N + 1;
        EXIT WHEN (Real (N*N) >= (3.2880242E15 / Frequency));
      END LOOP;
      N_prime := Integer'Min (8,N+2);
    END IF;

    --  Strom's fit to Karsas and Latter f-f Gaunt factors

    IF  (Frequency < 5.5E14)  THEN
      G3 := 1.0828 + 3.865E-6 * T + (7.564E-7 + (4.920E-10 - 2.482E-15 * T)
              * T + (5.326E-12 + (-3.904E-15 + 1.8790E-20 * T) * T) * W) * W;

    -- wavelength below 5450A
    ELSE
      G3 := 1.0;
    END IF;

    --  Kalkofen's integral for long wavelengths
    IF  (N_prime < N)  THEN
      D := (1.0 / Ehvkt - 1.0 + G3) / 2.0 / U;

    --  summation of contributions from N_prime to N relevant levels
    --  plus asymptotic relation for levels N_prime+1 and higher
    --  numerical constants from Allen, "Astrophysical Quantities", p. 89
    ELSE
      FOR  K  IN  N..N_prime  LOOP
        Square := Real (K*K);
        IF  (K = 1)  THEN
          B := Ecart;
        ELSE
          B := 1.0;
        END IF;
        D := B * Gaunt (Frequency, K) / Square * Exp (U / Square) / Real (K) + D;
      END LOOP;
      D := D + (Exp (U / Real ((N_prime+1)**2)) - 1.0 + G3) / 2.0 / U;
    END IF;
    RETURN  (1.0E26 / Frequency) * ( D / Exp (U)) * ((1.0 - Ehvkt) /
             Frequency * 2.815E29 / Frequency) / Ecart;
  END  Opac_Hydrogen;

  -------------------------------------------------------
  FUNCTION  Opac_H2p (Theta     : Real;
                      Frequency : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  H2+ opacity per neutral hydrogen atom and per H+, scaled by E26.
    --  Based on Bates, J.CHEM.PHYS. 19, 1122, 1951 -- MNRAS 112, 40, 1952
    --  Phil. Trans. R. London A 246, 215, 1953.
    --  see Gingerich, SAO SPECIAL REPORT 167, 21, 1964.
    --  12 February 1962 - 29 January 1964.

    N               : Integer  := 1;
    Ev              : Real     := 3.03979E-16 * Frequency;
    a0, a1, a2, a3  : Real;

    Elim : Const_Arr :=
       (1.0E36,   2.716,    2.437,   2.156,   1.898,    1.660,    1.445,
        1.254,    1.086,    0.9354,  0.8109,  0.7004,   0.6051,   0.5229,
        0.4518,   0.3904,   0.3374,  0.2915,  0.2517,   0.2171,   0.1871,
        0.1611,   0.1386,   0.1190,  0.1021,  0.08758,  0.07498,  0.06390,
        0.05434,  0.04629,  0.03950, 0.03357, 0.02840,  0.02397,  0.02029,
        0.01722,  0.01456,  0.01225, 0.01028, 0.008657, 0.007312, 0.006153,
        0.005149, 0.004300, 0.0    , 0.0      );

    E : Const_Arr :=
       (3.0,     2.852,   2.58,    2.294,   2.023,   1.774,   1.547,
        1.344,   1.165,   1.007,   0.8702,  0.7516,  0.6493,  0.5610,
        0.4848,  0.4189,  0.3620,  0.3128,  0.2702,  0.2332,  0.2011,
        0.1732,  0.1491,  0.1281,  0.1100,  0.09426, 0.0809,  0.06906,
        0.05874, 0.04994, 0.04265, 0.03635, 0.0308,  0.026,   0.02195,
        0.01864, 0.01581, 0.01332, 0.01118, 0.00938, 0.00793, 0.00669,
        0.00561, 0.00469, 0.00392, 0.0033  );

    Up : Const_Arr :=
       (85.0,        9.99465,     4.97842,     3.28472,     2.41452,
        1.87038,     1.48945,     1.20442,     0.98279,     0.80665,
        0.66493,     0.54997,     0.45618,     0.37932,     0.31606,
        0.26382,     0.22057,     0.18446,     0.15473,     0.12977,
        1.08890E-01, 9.14000E-02, 7.67600E-02, 6.44500E-02, 5.41200E-02,
        4.54000E-02, 3.81000E-02, 3.19500E-02, 2.67600E-02, 2.23700E-02,
        1.86900E-02, 1.56100E-02, 1.30200E-02, 1.08300E-02, 8.99000E-03,
        7.45000E-03, 6.15000E-03, 5.08000E-03, 4.16000E-03, 3.42000E-03,
        2.77000E-03, 2.21000E-03, 1.78000E-03, 1.45000E-03, 1.24000E-03,
        1.14000E-03  );

    Us : Const_Arr :=
       (-85.0,       -7.1426,     -2.3984,     -0.99032,    -0.39105,
        -0.09644,     0.05794,     0.13996,     0.18186,     0.20052,
         0.20525,     0.20167,     0.19309,     0.18167,     0.16871,
         0.15511,     0.14147,     0.12815,     0.11542,     0.10340,
         0.09216,     8.18000E-02, 7.22900E-02, 6.36700E-02, 5.58400E-02,
         4.88400E-02, 4.25700E-02, 3.69900E-02, 3.20700E-02, 2.77500E-02,
         2.39400E-02, 2.06100E-02, 1.77200E-02, 1.52200E-02, 1.30500E-02,
         1.11900E-02, 9.58000E-03, 8.21000E-03, 7.01000E-03, 6.00000E-03,
         5.11000E-03, 4.35000E-03, 3.72000E-03, 3.22000E-03, 2.86000E-03,
         2.63000E-03   );

    Fr : Const_Arr :=
       (0.0,         4.30272E-18, 1.51111E-17, 4.02893E-17, 8.89643E-17,
        1.70250E-16, 2.94529E-16, 4.77443E-16, 7.25449E-16, 1.06238E-15,
        1.50501E-15, 2.08046E-15, 2.82259E-15, 3.76256E-15, 4.93692E-15,
        6.38227E-15, 8.17038E-15, 1.02794E-14, 1.28018E-14, 1.57371E-14,
        1.91217E-14, 2.30875E-14, 2.75329E-14, 3.27526E-14, 3.85481E-14,
        4.52968E-14, 5.18592E-14, 5.99825E-14, 6.92092E-14, 7.94023E-14,
        9.01000E-14, 1.01710E-13, 1.14868E-13, 1.29969E-13, 1.46437E-13,
        1.63042E-13, 1.81440E-13, 2.02169E-13, 2.25126E-13, 2.49637E-13,
        2.73970E-13, 3.00895E-13, 3.30827E-13, 3.64140E-13, 3.99503E-13,
        4.34206E-13  );

  BEGIN
    LOOP
      EXIT WHEN (Ev >= Elim(N));
      N := N + 1;
    END LOOP;
    a0 := 0.0319273 / Theta;
    a1 := Para (X => E, Y => Fr, Arg_x => Ev, N => N);
    a2 := Para (X => E, Y => Us, Arg_x => Ev, N => N);
    a3 := Para (X => E, Y => Up, Arg_x => Ev, N => N);
    RETURN  ABS (a1 * (Exp (a2 / a0) - Exp (-a3 / a0)));
  END  Opac_H2p;

  -----------------------------------------------------------
  FUNCTION  Opac_H_minus (Theta     : Real;
                          Frequency : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  H- opacity per neutral hydrogen atom and unit electron presssure,
    --  scaled by E26. All opacity routines include stimulated emission.
    --  f-f based on John, MNRAS 128,93,1964.
    --  b-f based on Geltman, ApJ 136,935,1962.
    --  see Gingerich, SAO SPECIAL REPORT 167, 19-20, 1964.

    W  : Real  := cLight * 1.0E-3 / Frequency;
    X  : Real  := 16.419 - W;
    H  : Real  := 9.5210E-15 * Frequency;
    L  : Real;
  BEGIN

    --  f-f only
    IF  (X < 0.0)  THEN
      RETURN  0.0053666 + (-0.011493 + 0.027039 * Theta) * Theta +
              (-3.2062 + (11.924 - 5.939 * Theta) * Theta +
              (-0.40192 + (7.0355 - 0.34592 * Theta) * Theta) * W) * W / 1000.0;

    --  b-f plus f-f
    ELSE
      IF  (W > 14.2)  THEN
        L := (0.269818 + (0.220190 + (-0.0411288 + 0.00273236 * X) *
              X) * X) * X;
      ELSE
        L := 0.00680133 + (0.178708 + (0.16479 + (-0.0204842 + 5.95244E-4 *
              W) * W) * W) * W;
      END IF;
      RETURN  0.0053666 + (-0.011493 + 0.027039 * Theta) * Theta +
           (-3.2062 + (11.924 - 5.939 * Theta) * Theta +
           (-0.40192 + (7.0355 - 0.34592 * Theta) * Theta) * W) * W / 1000.0 +
            L * 0.4158 * Theta**2 * Sqrt (Theta) *
            Exp (1.737 * Theta) * (1.0 - Exp (-Theta * H));
    END IF;
  END  Opac_H_minus;

  ---------------------------------------------------------
  FUNCTION  Opac_Lyman (Theta     : Real;
                        Frequency : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Opacity per neutral hydrogen atom for resonance broadening of
    --  Lyman alpha line, based on calculations of Kenneth Sando, July 69.

    W : Real  := cLight / Frequency;
  BEGIN
    IF  ((W > 1950.0)  OR  (W < 1630.0))  THEN
      RETURN  0.0;
    ELSE

    --  Q.M. treatment for 1950A - 1630A .
      RETURN  4.03E-12 / Sqrt (5040.4 / Theta) *
                         Exp ((0.0471 + Theta * 35.0 / 5040.4) * (1620.0 - W));
    END IF;
  END  Opac_Lyman;

  --------------------------------------------------------------
  FUNCTION  Opac_Rayleigh_H (Frequency : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Rayleigh scattering per neutral hydrogen atom, scaled by E26.
    --  Dalgarno (see Gingerich, SAO SPECIAL REPORT 167, 20, 1964)

    W : Real  := (cLight / Real'Min (Frequency, 29.22E14))**2;
  BEGIN
    RETURN  (5.799E13 * (W**2) + 1.422E20 * W + 2.78E26) / W**4;
  END  Opac_Rayleigh_H;

  ---------------------------------------------------------------
  FUNCTION  Opac_Rayleigh_H2 (Frequency : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Rayleigh scattering per H2 molecule, scaled by E26 .
    --  Dalgarno & Williams (1962) - ApJ 136, p. 690 .

    W :  Real := (Frequency / cLight)**2;
  BEGIN
    RETURN  ((1.61E26 * W + 1.28E20) * W + 8.14E13) * W**2;
  END  Opac_Rayleigh_H2;

  ------------------------------------------------------------------
  FUNCTION  Opac_Quasi_Hydrogen (Theta     : Real;
                                 Frequency : Real)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --   Quasi_Hydrogen approximates the Quasi_Hydrogen molecular absorption.
    --   Coefficients given in the Doyle thesis, Harvard 1964.
    --   4 April 1968 version must be multiplied by H/cLight.
    --   To obtain opacity per neutral H Atom, multiply Quasi-Hydrogen by N(HI),
    --   number of HI atoms per cc..
    --   When multiplied by N(HI)**2, Quasi-Hydrogen gives opacity per unit volume.

    T          : Real                            := 5.0404 / Theta;
    Freq_limit : CONSTANT ARRAY (1..5)  OF Real :=
           (18.224E14, 15.989E14, 14.276E14, 12.757E14, 10.901E14);

    --   for     1645A     1875A      2100A      2350A      2750A

  BEGIN
    IF  (Frequency < Freq_limit(5))  THEN
      RETURN  0.0;
    ELSIF  (Frequency >= Freq_limit(1))  THEN                         ---  1540A
      RETURN  (-1.5527E3 + (2.0309E3 - 1.2391E2 * T) * T) * 1.0E-19;
    ELSIF  (Frequency >= Freq_limit(2))  THEN                         ---  1750A
      RETURN  (-2.2983E3 + (1.1117E3 - 40.873 * T) * T) * 1.0E-19;
    ELSIF  (Frequency >= Freq_limit(3))  THEN                         ---  2000A
      RETURN  (-462.84 + (140.64 + 11.809 * T) * T) * 1.0E-19;
    ELSIF  (Frequency >= Freq_limit(4))  THEN                         ---  2200A
      RETURN  (-24.813 + (-32.670 + 14.373 * T) * T) * 1.0E-19;
    ELSE                                                              ---  2500A
      RETURN  (79.057 + (-57.532 + 10.191 * T) * T) * 1.0E-19;
    END IF;
  END  Opac_Quasi_Hydrogen;

  ------------------------------------------------------
  FUNCTION  Chi_C (Temp   : Real;
                   Wave   : Real;
                   Zeta_0 : Real;
                   Bl     : Level_Arr)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Opacity per unit volume (cm-1) due to the photoionisation of
    --  neutral carbon (C I) at the wavelength W (in Angstroems).
    --  Zeta_0 = N0 / U0 = number of C I atoms per cm3 divided by the
    --  partition function of C I .
    --  Bl(L) = Departure coefficients from LTE of the population of level L .
    --  One takes into account 8 levels of C I .
    --  Numerical data -- see Vernazza et al., ApJS, 30, p. 6, Table 5 (1976) .

    W_lim : CONSTANT ARRAY (1..8)  OF Real :=
        ( 1100.0, 1239.0, 1444.0, 1745.0, 3257.0, 3437.0, 3733.0, 4907.0 );
    Al : CONSTANT ARRAY (1..8)  OF Real :=
         ( 40.4, 28.7, 33.6, 1.0,  0.2, 1.54,  16.0,  2.1 );
    Sl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 2.0, 1.5, 1.5, 3.0, 1.2, 1.2, 3.0, 1.5 );
    Gl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 9.0, 5.0, 1.0, 5.0, 9.0, 3.0, 15.0, 27.0 );
    Alp : CONSTANT ARRAY (1..8)  OF Real :=
         ( 28.2, 18.4, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    Slp : CONSTANT ARRAY (1..8)  OF Real :=
         ( 3.2, 2.5, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0 );
    Dw : CONSTANT ARRAY (1..8)  OF Real :=
         ( 0.0,         1.019884E-4, 2.165701E-4, 3.360250E-4,
           6.020599E-4, 6.181395E-4, 6.412098E-4, 7.053004E-4 );

    LM   : CONSTANT Integer := 8;
    Opac : Real             := 0.0;
    Sig  : Real;
  BEGIN
    FOR  L  IN  1..LM  LOOP
      IF  (Wave <= W_lim(L))  THEN
        Sig := Al(L) * (Wave / W_lim(L))**Sl(L);
        IF  (L <= 3)  THEN
          Sig := Sig - Alp(L) * (Wave / W_lim(L))**Slp(L);
        END IF;
        Opac := Opac + Sig * (1.0 - Exp (-C2 / (Wave * Temp)) / Bl(L)) *
                       (Zeta_0 * Gl(L) * Exp (-C2 * Dw(L) / Temp) * Bl(L));
      END IF;
    END LOOP;
    RETURN  Opac * 1.0E8;
  END  Chi_C;

  -------------------------------------------------------
  FUNCTION  Chi_Mg (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Opacity per unit volume (cm-1) due to the photoionisation of
    --  neutral magnesium (Mg I) at the wavelength W (in Angstroems).
    --  Zeta_0 = N0 / U0 = number of Mg I atoms per cm3 divided by the
    --  partition function of C I .
    --  Bl(L) = Departure coefficients from LTE of the population of level L .
    --  One takes into account 8 levels of Mg I .
    --  Numerical data -- see Vernazza et al., ApJS, 30, p. 6, Table 6 (1976) .

    W_lim : CONSTANT ARRAY (1..8)  OF Real :=
        ( 1621.0, 2513.0, 3756.0, 4884.0, 5504.0, 6549.0, 7236.0, 7292.0 );
    Al : CONSTANT ARRAY (1..8)  OF Real :=
         ( 1.1, 20.0, 16.0, 2.1, 0.43, 45.0, 25.0, 33.8 );
    Sl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 10.0, 2.7, 2.1, 2.6, 2.6, 2.7, 2.7, 2.8 );
    Gl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 1.0, 9.0, 3.0, 3.0, 1.0, 5.0, 9.0, 15.0 );
    Dw : CONSTANT ARRAY (1..8)  OF Real :=
         (0.0,         2.189724E-4, 3.506625E-4, 4.121529E-4,
          4.352171E-4, 4.642081E-4, 4.787052E-4, 4.797666E-4 );

    --   DW(L) = (1./WLIM(1)) - (1./WLIM(L))

    LM   : CONSTANT Integer := 8;
    Sig3 : Real             := 7.8 * (Wave / W_lim(3))**9.5;
    Opac : Real             := 0.0;
    Sig  : Real;
  BEGIN
    FOR  L  IN  1..LM  LOOP
      IF  (Wave <= W_lim(L))  THEN
        Sig := Al(L) * (Wave / W_lim(L))**Sl(L);
        IF  (L = 3)  THEN
          Sig := Sig - Sig3;
        END IF;
        Opac := Opac + Sig * (1.0 - Exp (-C2 / (Wave * Temp)) / Bl(L)) *
                      (Zeta_0 * Gl(L) * Exp (-C2 * Dw(L) / Temp) * Bl(L));
      END IF;
    END LOOP;
    RETURN  Opac * 1.0E8;
  END  Chi_Mg;

  -------------------------------------------------------
  FUNCTION  Chi_Al (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Opacity per unit volume (cm-1) due to the photoionisation of
    --  neutral aluminum (Al I) at the wavelength W (in Angstroems).
    --  Zeta_0 = N0 / U0 = number of C I atoms per cm3 divided by the
    --  partition function of C I .
    --  Bl(L) = Departure coefficients from LTE of the population of level L .
    --  One takes into account 8 levels of Al I .
    --  Numerical data -- see Vernazza et al., ApJS, 30, p. 6, Table 7 (1976) .

    W_lim : CONSTANT ARRAY (1..8)  OF Real :=
         (2076.0, 4360.0, 5205.0, 6311.0, 6525.0, 9442.0, 10698.0, 12495.0 );
    Al : CONSTANT ARRAY (1..8)  OF Real :=
         ( 65.0, 10.0, 10.0, 47.0, 14.5, 56.7, 50.0, 50.0 );
    Sl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 4.4, 2.0, 2.0, 1.83, 1.0, 1.9, 3.0, 3.0 );
    Gl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 6.0, 2.0, 6.0, 10.0, 6.0, 2.0, 10.0, 6.0 );
    Dw : CONSTANT ARRAY (1..8)  OF Real :=
         (0.0,         2.523378E-4, 2.895726E-4, 3.232421E-4,
          3.284389E-4, 3.757858E-4, 3.882202E-4, 4.016636E-4 );

    LM   : CONSTANT Integer := 8;
    Opac : Real             := 0.0;
    Sig  : Real;
  BEGIN
    FOR  L  IN  1..LM  LOOP
      IF  (Wave <= W_lim(L))  THEN
        Sig  := Al(L) * (Wave / W_lim(L))**Sl(L);
        Opac := Opac + Sig * (1.0 - Exp (-C2 / (Wave * Temp)) / Bl(L)) *
                      (Zeta_0 * Gl(L) * Exp (-C2 * Dw(L) / Temp) * Bl(L));
      END IF;
    END LOOP;
    RETURN  Opac * 1.0E8;
  END  Chi_al;

  -------------------------------------------------------
  FUNCTION  Chi_Si (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Opacity per unit volume (cm-1) due to the photoionisation of
    --  neutral silicon (Si I) at the wavelength W (in Angstroems).
    --  Zeta_0 = N0 / U0 = number of C I atoms per cm3 divided by the
    --  partition function of C I .
    --  Bl(L) = Departure coefficients from LTE of the population of level L .
    --  One takes into account 8 levels of Si I .
    --  Numerical data -- see Vernazza et al., ApJS, 30, p. 4, Table 1 (1976) .

    W_lim : CONSTANT ARRAY (1..8)  OF Real :=
         ( 1525.0, 1682.0, 1986.0, 3085.0, 3864.0, 4040.0, 4892.0, 5840.0 );
    Al : CONSTANT ARRAY (1..8)  OF Real :=
         ( 37.0, 35.0, 46.0, 15.0, 1.25, 4.09, 18.0, 14.1 );
    Sl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 5.0, 3.0, 0.5, 3.0, 2.0, 2.0, 3.0, 3.0 );
    Gl : CONSTANT ARRAY (1..8)  OF Real :=
         ( 9.0, 5.0, 1.0, 5.0, 9.0, 3.0, 15.0, 27.0 );
    Dw : CONSTANT ARRAY (1..8)  OF Real :=
         (0.0,         6.120738E-5, 1.522130E-4, 3.315886E-4,
          3.969385E-4, 4.082130E-4, 4.513223E-4, 4.845048E-4 );

    --  DW(L) = (1./WLIM(1)) - (1./WLIM(L))

    LM   : CONSTANT Integer  := 8;
    Opac : Real              := 0.0;
    Sig1 : Real              := Al(1);
    Sig  : Real;
  BEGIN
    IF  (Wave < 1350.0)  THEN
      Sig1 := Sig1 * (Wave / 1350.0)**SL(1);
    END IF;
    FOR  L  IN  1..LM  LOOP
      IF  (Wave <= W_lim(L))  THEN
        Sig := Al(L) * (Wave / W_lim(L))**Sl(L);
        IF  (L = 1)  THEN
          Sig := Sig1;
        END IF;
        Opac := Opac + Sig * (1.0 - Exp (-C2 / (Wave * Temp)) / Bl(L)) *
                      (Zeta_0 * Gl(L) * Exp (-C2 * Dw(L) / Temp) * Bl(L));
      END IF;
    END LOOP;
    RETURN  Opac * 1.0E8;
  END  Chi_Si;

  -------------------------------------------------------
  FUNCTION  Chi_Fe (Temp   : Real;
                    Wave   : Real;
                    Zeta_0 : Real;
                    Bl     : Level_Arr)  RETURN  Real  IS

    PRAGMA OPTIMIZE (Time);

    --  Opacity per unit volume (cm-1) due to the photoionisation of
    --  neutral Iron (Fe I) at the wavelength W (in Angstroems).
    --  Zeta_0 = N0 / U0 = number of C I atoms per cm3 divided by the
    --  partition function of C I .
    --  Bl(L) = Departure coefficients from LTE of the population of level L .
    --  One takes into account 2 levels of Fe I .
    --  Numerical data -- see Vernazza et al., ApJS, 30, p. 7, Table 8 (1976) .

    W_lim : CONSTANT ARRAY (1..2)  OF Real :=
         ( 1570.0, 1761.0 );
    Al : CONSTANT ARRAY (1..2)  OF Real :=
         ( 6.3, 5.04 );
    Sl : CONSTANT ARRAY (1..2)  OF Real :=
         ( 3.0, 3.0 );
    Gl : CONSTANT ARRAY (1..2)  OF Real :=
         ( 9.0, 11.0 );
    Dw : CONSTANT ARRAY (1..2)  OF Real :=
         (0.0, 6.908350E-5 );

    LM   : CONSTANT Integer  := 2;
    Opac : Real              := 0.0;
    Sig  : Real;
  BEGIN
    FOR  L  IN  1..LM  LOOP
      IF  (Wave <= W_lim(L))  THEN
        Sig  := Al(L) * (Wave / W_lim(L))**Sl(L);
        Opac := Opac + Sig * (1.0 - Exp (-C2 / (Wave * Temp)) / Bl(L)) *
                      (Zeta_0 * Gl(L) * Exp (-C2 * Dw(L) / Temp) * Bl(L));
      END IF;
    END LOOP;
    RETURN  Opac * 1.0E8;
  END  Chi_Fe;

END  Continuous_Opacity;

GENERIC
  TYPE  Real      IS DIGITS  <>;
  TYPE  Real_Arr  IS ARRAY (Integer RANGE <>)  OF  Real;
PACKAGE  Provide_partition_function  IS

  --    Provides the partition functions for the Nma degrees
  --    of ionisation of the element Iel as a function of
  --    temperature T (Theta = 5040.0 / T) according to the
  --    procedure of Traving, Baschek and Holweger, Abh. der
  --    Hamburger Sternwarte, Vol. VIII, Nr. 1, 1966 .
  --
  --    Written 1996 by Martin J. Stift
  --                    Institut fuer Astronomie
  --                    Tuerkenschanzstr. 17
  --                    A-1180 Wien
  --                    AUSTRIA
  --
  --                    e-mail:  stift@astro.univie.ac.at
  --
  --    This program is free software; you can redistribute it and/or modify
  --    it under the terms of the GNU General Public License as published by
  --    the Free Software Foundation; either version 2 of the License, or
  --    (at your option) any later version.
  --
  --    This program is distributed in the hope that it will be useful,
  --    but WITHOUT ANY WARRANTY; without even the implied warranty of
  --    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  --    GNU General Public License for more details.
  --
  --    You should have received a copy of the GNU General Public License
  --    along with this program; if not, write to the Free Software
  --    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  TYPE  Atom_type  IS RECORD
    Isyb  : Integer;  -- Identification
    Eps_g : Real;  -- Abundances, IAU Comm. 12, Grenoble, 1976 (AB(HE) = 0.1).
    Eps_k : Real;  -- Abundances according to Kurucz (1979)
    Mass  : Real;  -- Atomic masses, Allen, Astrophysical Quantities (1973).
  END RECORD;

  No_elements : CONSTANT Integer      :=  21;
  Nma_max     : CONSTANT Integer      :=   8;
  Q           : Real_Arr(1..Nma_max);
  Chi_ion     : Real_Arr(1..Nma_max);

  Atoms : CONSTANT ARRAY (1..No_elements)  OF  Atom_type :=
      ( ( 1,  12.0, 12.0,   1.008),   ( 2, 11.0, 11.05,  4.0026),
        ( 6,  8.67,  8.57, 12.0111),  ( 7,  7.96, 8.06, 14.0067),
        ( 8,  8.90,  8.83, 15.9994),  (10,  7.90, 7.45, 20.179),
        (11,  6.29,  6.24, 22.9898),  (12,  7.56, 7.54, 24.305),
        (13,  6.40,  6.40, 26.9815),  (14,  7.60, 7.55, 28.086),
        (15,  5.45,  5.43, 30.9738),  (16,  7.25, 7.21, 32.06),
        (18,  6.30,  6.65, 39.948),   (19,  5.14, 5.05, 39.102),
        (20,  6.32,  6.33, 40.08),    (22,  5.00, 4.74, 47.90),
        (24,  5.70,  5.70, 51.996),   (25,  5.40, 5.20, 54.938),
        (26,  7.50,  7.55, 55.847),   (27,  5.00, 4.50, 58.9332),
        (28,  6.30,  6.28, 58.71 )                                  );

  -----------------------------------------
  FUNCTION  Exp10 (X : Real)  RETURN  Real;

  --------------------------------------------
  FUNCTION  Qas (L    : Real;
                 N    : Integer;
                 DChi : Real;
                 Th   : Real   )  RETURN Real;

  --------------------------------------------------------------------------
  PROCEDURE  Partition_function
   (Iflag   : IN OUT  Integer;     -- if successful => 1
    Iel     : IN      Integer;     -- number of element in periodic table
    DChi    : IN      Real;        -- lowering of ionisation potential
    Theta   : IN      Real;        -- inverse temperature
    Q       : IN OUT  Real_arr;    -- partition functions
    Chi_ion : IN OUT  Real_arr;    -- ionisation potentials
    Nma     : IN OUT  Integer  );  -- no. of ionsation stages considered

END  Provide_partition_function;
WITH  Ada.Numerics.Generic_Elementary_Functions;
PACKAGE BODY  Provide_partition_function  IS

  --    Provides the partition functions for the Nma degrees
  --    of ionisation of the element Iel as a function of
  --    temperature T (Theta = 5040.0 / T) according to the
  --    procedure of Traving, Baschek and Holweger, Abh. der
  --    Hamburger Sternwarte, Vol. VIII, Nr. 1, 1966 .
  --
  --    Written 1996 by Martin J. Stift
  --                    Institut fuer Astronomie
  --                    Tuerkenschanzstr. 17
  --                    A-1180 Wien
  --                    AUSTRIA
  --
  --                    e-mail:  stift@astro.univie.ac.at
  --
  --    This program is free software; you can redistribute it and/or modify
  --    it under the terms of the GNU General Public License as published by
  --    the Free Software Foundation; either version 2 of the License, or
  --    (at your option) any later version.
  --
  --    This program is distributed in the hope that it will be useful,
  --    but WITHOUT ANY WARRANTY; without even the implied warranty of
  --    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  --    GNU General Public License for more details.
  --
  --    You should have received a copy of the GNU General Public License
  --    along with this program; if not, write to the Free Software
  --    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  PACKAGE  Math  IS NEW  Ada.Numerics.Generic_Elementary_Functions (Real);
  USE      Math;

  --------------------------------------------
  FUNCTION  Exp10 (X : Real)  RETURN  Real  IS
    Ln10 : CONSTANT Real  := Log (10.0);
  BEGIN
    RETURN  Exp (x * Ln10);
  END;

  ----------------------------------------------
  FUNCTION  Qas (L    : Real;
                 N    : Integer;
                 DChi : Real;
                 Th   : Real   )  RETURN Real IS

    --  effective nuclear charge H

    cRydberg : CONSTANT Real := 13.595;         -- eV
    H        : CONSTANT Real := Sqrt (Real (N) * cRydberg / DChi);
    Effk     : CONSTANT Real := Th * 31.321 * Real (N**2);
  BEGIN
    RETURN  (H * (H + 1.0) * (H + 0.5) - L * (L - 1.0) * (L - 0.5))
             / 3.0 - Effk * (H - L + 1.0)
             * (1.0 - (Effk / 2.0) / H / (L - 1.0));
  END  Qas;

  ---------------------------------
  PROCEDURE  Partition_function
   (Iflag   : IN OUT  Integer;
    Iel     : IN      Integer;
    DChi    : IN      Real;
    Theta   : IN      Real;
    Q       : IN OUT  Real_arr;
    Chi_ion : IN OUT  Real_arr;
    Nma     : IN OUT  Integer )  IS

    TYPE  Part_data_type  IS ARRAY (1..175)  OF  Real;

    G : ARRAY (1..No_elements)  OF  Part_data_type :=

      ( 1 => (  2.000,   1.000,  2.000,   2.000,  2.000, 13.595, 11.000,
               10.853,  20.498, 13.342, 747.502,  1.000,  1.000,  0.000,
                0.000, 999.999,  1.000,
               OTHERS => 0.0                                                  ),

        2 => (  2.000,  1.000,  1.000,   2.000,  4.000, 24.580,   8.000,
               21.170, 28.170, 24.125, 527.830,  1.000,  2.000,   2.000,
                2.000, 54.403, 12.000,  43.708, 22.281, 53.542, 987.719,
               OTHERS => 0.0                                                  ),

        3 => (  5.000,   1.000,   1.000,  4.000, 12.000,  11.256,   6.000,
                0.004,   8.016,   1.359,  5.883,  6.454,  33.752,  10.376,
              595.343,   2.000,   2.000,  3.000,  2.000,  24.376,   6.000,
                0.008,   4.000,  16.546, 17.084, 21.614,  82.915,   3.000,
               18.000,  30.868,   5.000,  5.688, 15.981,  15.801,  48.204,
               26.269, 435.809,   2.000,  1.000,  3.000,   4.000,  47.871,
                6.100,   6.691,  10.028, 25.034, 15.757,  40.975, 186.211,
                3.000,  12.000,  55.873,  5.000, 17.604,  15.413,  36.180,
               55.956,  47.133, 243.631,  1.000,  2.000,   3.000,   2.000,
               64.476,   6.000,   8.005,  6.006, 40.804,  23.576,  54.492,
               76.418,   1.000,   1.000,  2.000,  4.000, 391.986,   4.000,
              303.772,  15.799, 354.208, 36.200,
              OTHERS => 0.0                                                   ),

        4 => (  5.000,   2.000,  4.000,  3.000,  18.000,  14.529,   6.100,
                2.554,  14.050,  9.169, 30.801,  13.651, 883.144,   3.000,
               10.000,  16.428,  4.000, 12.353,  10.000,  13.784,  16.000,
               14.874,  64.000,  2.000,  1.000,   4.000,  12.000,  29.593,
                5.000,   0.014,  8.046,  2.131,   6.267,  15.745,  17.870,
               24.949, 282.808,  3.000, 24.000,  36.693,   3.900,   6.376,
                7.375,  14.246, 33.139, 29.465, 215.483,   3.000,   2.000,
                3.000,   2.000, 47.426,  6.000,   0.022,   4.000,  31.259,
               19.353,  41.428, 80.646,  4.000,  18.000,  55.765,   5.000,
                7.212,  13.100, 15.228, 19.642,  34.387,  94.303,  46.708,
              370.954,   2.000,  6.000, 63.626,   4.000,  46.475,  16.000,
               49.468,  38.000,  2.000,  1.000,   3.000,   4.000,  77.450,
                6.000,   8.693, 10.329, 37.650,  14.502,  65.479, 187.162,
                2.000,  12.000, 87.445,  6.300,  61.155, 108.161,  79.196,
              191.838,   1.000,  2.000,  3.000,   2.000,  97.863,   6.000,
                9.999,   6.004, 60.991, 23.561,  82.262,  76.434,
              OTHERS => 0.0                                                   ),

        5 => (  5.000,   3.000,   5.000,   3.000,   8.000,  13.614,   8.000,
                0.022,   4.003,   2.019,   5.366,   9.812,  36.285,   2.000,
               20.000,  16.938,   6.000,  13.804, 131.022,  16.061, 868.978,
                2.000,  12.000,  18.630,   3.400,  14.293,  14.853,  16.114,
               93.147,   4.000,   4.000,   4.000,  18.000,  35.108,   6.000,
                3.472,  12.784,   7.437,   5.683,  22.579,  98.092,  32.035,
              829.440,   2.000,  10.000,  37.621,   5.000,  27.774,  50.988,
               33.678, 199.012,   3.000,   2.000,  40.461,   3.900,  28.118,
                2.000,  31.019,   6.000,  34.204,  10.000,   3.000,  10.000,
               42.584,   3.900,  30.892,  10.000,  33.189,  30.000,  36.181,
               50.000,   3.000,   1.000,   4.000,  12.000,  54.886,   6.000,
                0.032,   8.070,   2.760,   5.714,  35.328,  84.116,  48.277,
              529.093,   4.000,  24.000,  63.733,   4.900,   7.662,   5.661,
               16.786,  28.936,  42.657, 111.362,  54.522, 494.041,   2.000,
               20.000,  70.556,   4.000,  50.204,  45.525,  56.044, 134.475,
                5.000,   2.000,   3.000,   2.000,  77.394,   5.900,   0.048,
                4.000,  50.089,  21.294,  66.604,  78.706,   4.000,  18.000,
               87.609,   5.000,   8.954,  12.829,  18.031,  16.273,  57.755,
              123.658,  72.594, 327.240,   2.000,   6.000,  97.077,   4.900,
               68.388,  48.788,  82.397, 102.212,   2.000,  18.000, 103.911,
                4.000,  31.960,  20.006,  76.876, 161.990,   2.000,  10.000,
              106.116,   4.000,  75.686,  28.418,  80.388,  61.582,   2.000,
                1.000,   3.000,   4.000, 113.873,   6.000,  10.747,  10.556,
               52.323,  13.295,  94.976, 188.139,   3.000,  12.000, 125.863,
                6.000,  27.405,  14.656,  86.350, 129.492, 109.917, 470.851,
              OTHERS => 0.0                                                   ),

        6 => (  5.000,   2.000,   1.000,   2.000,   8.000,  21.559,   6.000,
               17.796,  34.508,  20.730, 365.492,   2.000,   4.000,  21.656,
                6.000,  17.879,  16.577,  20.855, 183.423,   2.000,   4.000,
                3.000,  18.000,  41.071,   5.000,   0.097,   2.001,  29.878,
               89.561,  37.221, 380.438,   2.000,  10.000,  44.274,   4.000,
               31.913,  26.447,  37.551,  63.553,   3.000,   5.000,   4.000,
                8.000,  63.729,   3.900,   0.092,   4.034,   3.424,   5.616,
               24.806,  11.518,  46.616,  72.827,   2.000,  20.000,  68.806,
                4.000,  45.643,  48.568,  54.147, 131.431,   2.000,  12.000,
               71.434,   4.000,  48.359,  31.171,  57.420,  76.829,   2.000,
                4.000,   4.000,  18.000,  97.162,   5.000,   5.453,  14.048,
               18.560,  13.308,  46.583,  52.790,  80.101, 467.849,   2.000,
               10.000, 100.917,   5.000,  70.337,  54.220,  85.789, 195.780,
                1.000,   1.000,   4.000,  12.000, 126.423,   4.000,   0.135,
                8.381,   5.497,   7.972,  26.121,  39.332,  86.665, 154.283,
               OTHERS => 0.0                                                  ),

        7 => (  5.000,   1.000,   2.000,   2.000,   2.000,   5.138,   7.000,
                2.400,  11.635,   4.552, 158.359,   2.000,   1.000,   2.000,
                8.000,  47.290,   4.000,  34.367,  21.045,  40.566,  50.955,
                2.000,   4.000,  47.459,   4.000,  34.676,  10.139,  40.764,
               25.861,   2.000,   4.000,   3.000,  18.000,  71.647,   4.000,
                0.170,   2.002,  44.554,  38.057,  57.142, 137.940,   2.000,
               10.000,  75.504,   4.000,  51.689,  28.311,  60.576,  61.689,
                3.000,   5.000,   4.000,   8.000,  98.880,   5.000,   0.152,
                4.033,   4.260,   5.856,  36.635,  18.179,  83.254, 208.914,
                2.000,  20.000, 104.778,   5.000,  72.561,  93.689,  89.475,
              406.309,   2.000,  12.000, 107.864,   5.000,  75.839,  60.428,
               92.582, 239.572,   5.000,   4.000,   4.000,  18.000, 138.597,
                5.000,   5.930,  10.317,   9.127,   5.697,  94.219, 121.073,
              115.690, 340.913,   2.000,  10.000, 142.980,   5.000,  97.783,
               63.792,  12.110, 186.208,   2.000,   2.000, 147.803,   4.000,
               92.769,   2.000, 108.400,  16.000,   3.000,  10.000, 151.427,
                4.000,  28.372,  17.927,  49.002,  18.704, 111.473,  91.369,
                2.000,  30.000, 163.906,   4.000, 110.226,  36.852, 125.690,
              233.147,
              OTHERS => 0.0                                                   ),

        8 => (  5.000,   1.000,   1.000,   3.000,   4.000,   7.644,   7.000,
                2.805,  10.744,   6.777, 291.506,   9.254,  53.749,   1.000,
                2.000,   3.000,   2.000,  15.031,   7.000,   4.459,   6.227,
                9.789,  31.129,  13.137, 132.644,   2.000,   1.000,   2.000,
                8.000,  80.117,   5.000,  57.413,  40.438,  71.252, 159.562,
                2.000,   4.000,  80.393,   5.000,  58.010,  20.384,  71.660,
               79.615,   2.000,   4.000,   3.000,  18.000, 109.294,   5.000,
                0.276,   2.001,  74.440, 106.898,  94.447, 343.101,   2.000,
               10.000, 113.799,   5.000,  54.472,  10.133,  95.858, 237.858,
                3.000,   5.000,   4.000,   8.000, 141.231,   5.000,   0.251,
                4.110,   5.370,   6.154,  51.461,  22.364, 121.819, 275.334,
                2.000,  20.000, 147.944,   5.000, 101.026, 117.038, 124.652,
              382.962,   2.000,  12.000, 151.493,   5.000, 104.636,  71.423,
              128.118, 228.576,
              OTHERS => 0.0                                                   ),

        9 => (  5.000,   2.000,   2.000,   3.000,   2.000,   5.984,   7.000,
                0.014,   4.001,   3.841,  11.780,   5.420, 142.218,   2.000,
               18.000,  10.634,   4.000,   3.727,  13.658,   8.833,  96.337,
                2.000,   1.000,   3.000,   4.000,  18.823,   7.000,   4.749,
               10.081,  11.902,  49.584,  16.719, 285.334,   2.000,  12.000,
               25.496,   4.000,  11.310,  14.687,  18.268,  59.312,   1.000,
                2.000,   3.000,   2.000,  28.441,   7.000,   6.751,   6.328,
               16.681,  29.509,  24.151, 134.163,   2.000,   1.000,   2.000,
                8.000, 119.957,   5.000,  83.551,  46.316, 104.787, 153.683,
                2.000,   4.000, 120.383,   5.000,  84.293,  22.990, 105.171,
               77.010,   1.000,   4.000,   4.000,  18.000, 153.772,   5.000,
                0.426,   2.000,  46.012,   2.360, 107.928, 151.214, 131.729,
              392.425,
              OTHERS => 0.0                                                   ),

       10 => (  5.000,   1.000,   1.000,   5.000,  12.000,   8.149,   6.100,
                0.020,   7.966,   0.752,   4.676,   1.614,   1.351,   5.831,
              123.227,   7.431, 443.780,   2.000,   2.000,   4.000,   2.000,
               16.339,   5.900,   0.036,   4.000,   8.795,   7.419,  11.208,
               24.175,  13.835,  60.406,   4.000,  18.000,  22.894,   5.000,
                5.418,  14.469,   7.825,  11.972,  14.440,  26.506,  19.412,
              269.052,   2.000,   1.000,   4.000,   4.000,  33.459,   5.000,
                6.572,   9.179,  11.449,   4.877,  18.424,  29.144,  25.457,
               52.800,   3.000,  12.000,  42.333,   5.000,  15.682,  13.267,
               27.010,  36.042,  34.599, 180.691,   1.000,   2.000,   3.000,
                2.000,  45.130,   7.000,   9.042,   6.484,  24.101,  27.685,
               37.445, 135.830,   2.000,   1.000,   2.000,   8.000, 166.725,
                4.900, 115.608,  36.982, 142.257, 163.018,   2.000,   4.000,
              167.357,   4.900, 118.377,  19.116, 143.084,  80.884,
              OTHERS => 0.0                                                   ),

       11 => (  5.000,  2.000,   4.000,  3.000,  18.000,  10.474,   5.000,
                1.514, 13.521,   5.575, 22.213,   9.247, 353.258,   2.000,
               10.000, 11.585,   5.000,  8.076,  10.000,  10.735, 150.000,
                1.000,  1.000,   4.000, 12.000,  19.720,   5.000,   0.043,
                8.024,  1.212,   5.808,  8.545,  51.754,  15.525, 252.400,
                1.000,  2.000,   4.000,  2.000,  30.156,   7.000,   0.074,
                4.002,  7.674,  20.798, 16.639,  62.419,  25.118, 200.779,
                1.000,  1.000,   3.000,  4.000,  51.354,   8.600,   8.992,
               11.741, 24.473,  63.512, 40.704, 179.742,   1.000,   2.000,
                3.000,  2.000,  65.007,  8.000,  11.464,   6.883,  33.732,
               32.778, 55.455, 228.337,
               OTHERS => 0.0                                                  ),

       12 => (  5.000,   3.000,   5.000,   4.000,  8.000,  10.357,   6.000,
                0.053,   3.961,   1.121,   5.078,  5.812,  15.094,   9.425,
              362.859,   2.000,  20.000,  12.200,  5.000,   8.936,  51.599,
               11.277, 268.400,   2.000,  12.000, 13.401,   5.000,   9.600,
               12.000,  12.551, 276.000,   2.000,  4.000,   4.000,  18.000,
               23.405,   5.000,   1.892,  11.438,  3.646,   5.513,  13.550,
              141.001,  19.376, 254.048,   2.000, 10.000,  24.807,   5.000,
               16.253,  33.052,  21.062, 126.948,  1.000,   1.000,   5.000,
               12.000,  35.047,   3.500,   0.043,  4.071,   0.123,   4.064,
                1.590,   5.724,  13.712, 144.638, 22.050, 106.491,   2.000,
                2.000,   4.000,   2.000,  47.292,  5.000,   0.118,   4.001,
                9.545,  19.281,  18.179,  27.599, 31.441,  35.118,   2.000,
               18.000,  57.681,  14.400,  30.664, 94.745,  56.150, 283.249,
                2.000,   1.000,   3.000,   4.000, 72.474,   5.000,  10.704,
               10.547,  27.075,  28.714,  50.599, 65.738,   1.000,  12.000,
               85.701,   4.000,  43.034,  24.000,
              OTHERS => 0.0                                                   ),

       13 => (  5.000,  2.000,   1.000,   2.000,   8.000,  15.755,  6.000,
               12.638, 43.662,  14.958, 324.337,   2.000,   4.000, 15.933,
                6.000, 12.833,  20.830,  15.139, 163.170,   2.000,  4.000,
                3.000, 18.000,  27.619,   5.100,   0.178,   2.003, 17.522,
              137.451, 23.584, 258.544,   2.000,  10.000,  29.355,  5.000,
               20.464, 62.813,  25.150, 149.187,   3.000,   5.000,  4.000,
                8.000, 40.899,   5.000,   0.151,   4.049,   1.561, 14.447,
               17.399, 46.823,  30.871, 124.665,   2.000,  20.000, 42.407,
                5.000, 24.684, 151.983,  33.978, 268.016,   2.000, 12.000,
               45.234,  5.000,  27.091, 101.130,  36.481, 150.869,  1.000,
                4.000,  4.000,  18.000,  59.793,   5.000,   2.810, 13.372,
                8.877,  8.653,  24.351,  60.461,  44.489, 285.507,  1.000,
                1.000,  4.000,  12.000,  75.002,   4.000,   0.144,  6.765,
                1.160,  4.768,  10.210,  12.863,  27.178,  54.526,
              OTHERS => 0.0                                                   ),

       14 => (  5.000,   1.000,   2.000,   3.000,  2.000,  4.339,   7.000,
                1.871,  12.978,   3.713, 148.667, 18.172,  6.349,   2.000,
                1.000,   2.000,   8.000,  31.810,  5.000, 21.185,  66.344,
               27.705, 101.655,   3.000,   4.000, 32.079,  5.000,   2.059,
                4.000,  23.709,  13.446,  28.542, 46.553,  3.000,   4.000,
                3.000,  18.000,  45.738,   6.000,  0.273,  2.017,  26.709,
              116.477,  39.640, 713.496,   2.000, 10.000, 47.768,   6.000,
               31.220,  63.591,  41.865, 396.408,  3.000,  2.000,  50.515,
                5.000,  29.955,   2.000,  37.557, 10.000, 42.862,  30.000,
                3.000,   5.000,   4.000,   8.000, 60.897,  6.000,   0.228,
                4.070,   2.274,   5.779,  21.703, 52.679, 50.191, 327.454,
                2.000,  20.000,  63.890,   5.000, 32.145, 62.860,  49.262,
              357.133,   2.000,  12.000,  65.849,  5.000, 34.155,  55.934,
               51.718, 196.065,   2.000,   4.000,  4.000, 18.000,  82.799,
                3.600,   3.043,  10.927,   5.479,  5.540, 20.547,  43.276,
               30.680,  76.256,   2.000,  10.000, 85.150,  4.000,  36.275,
               42.000,  47.345,  18.000,
              OTHERS => 0.0                                                   ),

       15 => (  4.000,  2.000,  1.000,  3.000,   4.000,   6.111,   5.900,
                2.050, 18.237,  3.349, 27.501,   5.321, 149.262,   2.000,
               20.000,  7.808,  6.000,  4.873,  94.524,   7.017, 705.471,
                1.000,  2.000,  3.000,  2.000,  11.868,   7.000,   1.769,
               11.871,  5.109, 14.071,  9.524, 106.055,   2.000,   1.000,
                2.000,  8.000, 51.207,  5.000,  27.271,  57.241,  41.561,
              110.757,  2.000,  4.000, 51.596,   4.900,  29.172,  29.812,
               42.140, 54.187,  2.000,  4.000,   3.000,  18.000,  67.181,
                5.000,  0.394,  2.018, 28.930,  97.578,  52.618, 282.394,
                2.000, 10.000, 69.536,  4.300,  38.593, 209.187,  49.646,
              252.813,
              OTHERS => 0.0                                                   ),

       16 => (  5.000,  3.000,   5.000,  6.000,  56.000,   6.818,  5.000,
                0.021,  7.089,   0.048,  8.919,   1.029,  17.563,  2.183,
              206.683,  4.109, 438.573,  5.785, 654.172,   4.000, 56.000,
                6.953,  4.700,   0.846, 38.046,   1.792,  69.627,  3.836,
              364.284,  5.787, 832.041,  3.000,  28.000,   7.411,  5.000,
                2.561, 98.856,   4.869, 57.993,   6.340, 442.150,  2.000,
                4.000,  6.000,  42.000, 13.635,   5.000,   0.023, 19.784,
                0.124, 32.064,   0.774, 37.089,   1.810, 110.668,  4.980,
              288.495,  9.585, 521.884,  3.000,  10.000,  14.685,  5.000,
                1.082, 10.000,   4.928, 34.000,  11.279, 120.000,  1.000,
                5.000,  5.000,  20.000, 28.137,   5.000,   0.041, 16.169,
                1.375, 22.355,   4.768, 24.165,  10.985,  83.513, 19.769,
              222.796,  1.000,   4.000,  4.000,   2.000,  43.236,  6.000,
                0.048,  6.002,  11.577,  4.618,  24.531,  25.264, 36.489,
               52.116,  1.000,   1.000,  2.000,  12.000, 100.083,  4.800,
               54.436, 12.000,  75.373,  8.000,
               OTHERS => 0.0                                                  ),

       17 => (  6.000,   3.000,  7.000,   3.000,  12.000,   6.763,   5.000,
                0.993,  30.184,  3.070,  79.285,   5.673, 149.529,   3.000,
               60.000,   8.285,  5.000,   3.339, 215.370,   4.801, 119.197,
                7.198, 741.432,  1.000,  40.000,   9.221,   5.000,   2.829,
              184.995,   2.000,  6.000,   4.000,  50.000,  16.493,   5.000,
                1.645,  46.619,  3.727, 160.136,   7.181, 488.045,  12.299,
              657.193,   4.000, 18.000,  18.662,   5.000,   2.902,  47.174,
                4.273, 267.027,  8.569, 441.132,  14.912, 150.665,   1.000,
                1.000,   4.000, 56.000,  30.950,   5.000,   0.047,  24.377,
                2.556, 122.836,  9.441, 285.509,  21.198, 794.165,   1.000,
                4.000,   4.000, 42.000,  49.580,   5.000,   0.078,  24.230,
                2.242,  75.026, 15.638, 172.945,  32.725, 543.651,   1.000,
                5.000,   4.000, 20.000,  73.093,   5.200,   0.103,  15.982,
                2.146,  17.680, 26.153,  95.200,  49.381, 225.095,   1.000,
                4.000,   3.000,  2.000,  90.595,   4.900,   0.119,   6.001,
               32.711,   6.873, 58.117,  15.125,
              OTHERS => 0.0                                                   ),

       18 => (  6.000,   3.000,  6.000,   3.000,  14.000,   7.432,   6.000,
                2.527,  53.911,  4.204,  81.393,   6.602, 546.694,   2.000,
               10.000,   8.606,  5.000,   4.155, 144.189,   7.321, 407.803,
                2.000,  50.000,  9.240,   6.000,   2.285,  45.618,   5.631,
              298.442,   2.000,  7.000,   4.000,  12.000,  15.636,   7.000,
                1.496,  22.638,  3.839,  93.842,   7.751, 183.937,  13.484,
              907.576,   4.000, 72.000,  18.963,   5.000,   3.681, 137.041,
                6.054, 168.678,  9.934, 329.029,  14.936, 773.251,   1.000,
                6.000,   4.000, 50.000,  33.690,   5.000,   3.531,  70.192,
                6.967,  72.337, 15.222, 213.951,  25.069, 539.516,   1.000,
                1.000,   4.000, 56.000,  53.001,   5.000,   0.071,  24.237,
                2.896,  93.541, 20.725, 456.617,  37.383, 506.548,   1.000,
                4.000,   4.000, 42.000,  76.006,   5.000,   0.126,  24.769,
                2.660,  66.990, 28.528, 264.185,  53.413, 484.016,   1.000,
                5.000,   4.000, 20.000, 109.002,   3.700,   0.117,  12.314,
                0.350,   4.305,  2.539,  22.242,  40.301,  60.120,
              OTHERS => 0.0                                                   ),

       19 => (  5.000,   3.000,   9.000,   4.000,  60.000,   7.896,   5.000,
                0.066,  14.410,   0.339,   2.705,   2.897, 421.661,   6.585,
              940.148,   4.000,  56.000,   8.195,   5.000,   0.923,  36.219,
                1.679,  22.888,   4.620, 239.600,   7.053, 825.292,   3.000,
               40.000,   8.927,   5.000,   4.249, 110.024,   5.875, 992.304,
                7.781, 640.671,   2.000,  10.000,   4.000,  50.000,  16.178,
                5.000,   0.062,  17.049,   0.283,  32.378,   1.504,  34.318,
                5.430, 420.963,   3.000,  18.000,  18.662,   5.000,   2.792,
              154.006,   7.627, 462.112,  13.623, 329.862,   2.000,   9.000,
                4.000,  12.000,  30.640,   6.000,   0.077,  15.791,   3.723,
               47.119,  12.137, 279.929,  23.700, 692.100,   4.000,  72.000,
               34.607,   5.000,   2.688,  91.021,   7.595, 206.308,  15.444,
              706.993,  25.587, 836.669,   1.000,   6.000,   4.000,  50.000,
               56.001,   3.600,   3.982,  40.079,   4.677,  27.696,   6.453,
               28.224,  23.561,  18.000,   1.000,   1.000,   4.000,  56.000,
               79.001,   3.800,   0.102,  24.090,   3.354,  89.634,  22.954,
               51.576,  33.796, 241.698,
              OTHERS => 0.0                                                   ),

       20 => (  4.000,   4.000,  10.000,   5.000,  42.000,   7.863,   5.000,
                0.112,  11.912,   0.341,  20.442,   0.809,  28.386,   3.808,
              132.504,   6.723, 600.746,   3.000,  70.000,   8.378,   5.000,
                2.057,  33.309,   3.484, 237.433,   7.210, 977.250,   3.000,
               42.000,   9.160,   5.000,   2.405,  55.540,   5.133, 318.817,
                8.097, 619.637,   3.000,  18.000,   9.519,   5.000,   2.084,
               32.690,   5.291,  83.869,   8.426, 107.438,   2.000,   9.000,
                5.000,  56.000,  17.052,   5.000,   0.135,  11.259,   0.517,
               38.224,   1.606,  22.996,   6.772, 261.349,  12.622, 637.148,
                4.000,  24.000,  18.958,   5.000,   2.512,  23.023,   4.348,
               41.660,   8.253, 264.646,  15.377, 181.670,   1.000,  10.000,
                5.000,  50.000,  33.491,   5.000,   0.132,  16.036,   0.863,
                7.863,   3.086,  70.316,  11.789, 423.351,  23.263, 742.355,
                1.000,  20.000,   0.000,  12.000,  53.001,   3.000,
              OTHERS => 0.0                                                   ),

       21 => (  4.000,  2.000,   9.000,   6.000,  20.000,   7.633,  5.400,
                0.026,  7.127,   0.137,  12.449,   0.315,  11.995,  1.778,
               10.055,  4.029, 114.166,   6.621, 391.206,   3.000, 56.000,
                8.793,  5.000,   2.249,  26.391,   4.042, 213.808,  7.621,
              938.793,  2.000,   6.000,   5.000,  42.000,  18.147,  9.000,
                0.191,  4.142,   1.235,  37.378,   3.358,  25.971,  8.429,
              333.340, 17.096, 311.163,   3.000,  18.000,  20.233,  5.000,
                3.472, 33.103,   9.065, 184.185,  16.556, 136.707,  1.000,
                9.000,  5.000,  56.000,  35.161,   5.000,   0.194, 11.191,
                1.305,  5.417,   5.813,  53.679,  14.172, 460.678, 26.169,
              380.006,  1.000,  28.000,   0.000,  50.000,  56.025,  3.000,
              OTHERS => 0.0                                                  ));

    --  Calculate the partition functions for the Nma dgrees of ionisation of
    --  the element Iel as a function of temperature T (Theta = 5040.0 / T)
    --  according to the procedure of Traving, Baschek and Holweger,
    --  Abhandlungen der Hamburger Sternwarte, Vol. VIII, Nr. 1, 1966 .

    Ic, It, K         : Integer;
    Igg, Igp, G0, Gpr : Real;
    Sum, X, Y         : Real;

    Chi  : ARRAY (1..8)  OF  Real;
    GG   : ARRAY (1..8)  OF  Real;

  BEGIN
    Iflag := -IEL;
    FOR  M  IN  1..No_elements  LOOP
      Ic := 0;
      IF  (Iel = Atoms(M).Isyb)  THEN
        Nma := Integer (G(M)(Ic+1));
        IF  (Nma > 3)  THEN
          Nma := 3;
        END IF;
        Ic  := Ic + 1;
        FOR  N  IN  1..Nma  LOOP
          It   := Integer (G(M)(Ic+1));
          Igg  := G(M)(Ic+2);
          Ic   := Ic + 2;
          G0   := Igg;
          Q(N) := G0;
          FOR  I  IN  1..It  LOOP
            K    := Integer (G(M)(Ic+1));
            Igp  := G(M)(Ic+2);
            X    := G(M)(Ic+3);
            Y    := G(M)(Ic+4);
            Ic   := Ic + 4;
            Gpr  := Igp;
            IF  (I = 1)  THEN
              Chi_ion(N) := X;
            END IF;
            IF  (K /= 0)  THEN
              FOR  L  IN  1..K  LOOP
                Chi(L) := G(M)(Ic+1);
                GG(L)  := G(M)(Ic+2);
                Ic     := Ic + 2;
              END LOOP;
              Sum := 0.0;
              FOR  J  IN  1..K  LOOP
                Sum := Sum + GG(J) * Exp10 (-Chi(J) * Theta);
              END LOOP;
              Q(N) := Q(N) + Sum + Gpr * Exp10 (-Theta * X) *
                      Qas (L => Y, N => N, DChi => DChi, TH => -Theta);
            END IF;
          END LOOP;
        END LOOP;
        Iflag := 1;
      END IF;
    END LOOP;

  END  Partition_function;

END  Provide_partition_function;
