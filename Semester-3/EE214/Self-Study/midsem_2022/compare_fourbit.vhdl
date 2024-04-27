library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity compare_fourbit is
port (X, Y : in std_logic_vector(3 downto 0); 
		Z, N, V, L, G : out std_logic);
end entity compare_fourbit;

architecture Struct of compare_fourbit is
--TO BE CHECKED!
signal S3, S2, S1, S0, car, NY3, NX3, NS3, V1, V2, Zzbar, Zz, Nn, Vv, VxorredN, L1, Ll, nL, Gg: std_logic := '0';

component Fourbit_AdderSubtractor is
port (A3, A2, A1, A0, B3, B2, B1, B0, M: in std_logic; Cout, S3, S2, S1, S0 : out std_logic);
end component;

begin
  --component instances
  Subtract: Fourbit_AdderSubtractor PORT MAP (A3 => X(3), A2 => X(2), A1 => X(1), A0 => X(0), 
													B3 => Y(3), B2 => Y(2), B1 => Y(1), B0 => Y(0), M => '1', 
													Cout => car, S3 => S3, S2 => S2, S1 => S1, S0 => S0);
	Y3Bar : INVERTER PORT MAP (A => Y(3), Y => NY3);
	X3Bar : INVERTER PORT MAP (A => X(3), Y => NX3);
	S3Bar : INVERTER PORT MAP (A => S3, Y => NS3);
	checkEq1 : OR_4 PORT MAP (A => S3, B => S2, C => S1, D => S0, Y => Zzbar);
	checkEq2: INVERTER PORT MAP (A => Zzbar, Y => Zz);
	checkNEG : AND_2 PORT MAP (A => S3, B => S3, Y => Nn);
	checkOvF1 : AND_3 PORT MAP (A => S3, B => NX3, C => Y(3), Y => V1);	 --This, alone gives 28 errors (OvF > 7 eliminated)
	checkOvF2 : AND_3 PORT MAP (A => X(3), B => NS3, C => NY3, Y => V2);  --This, alone gives 36 errors (OvF < -8 eliminated)
	checkOvF3 : OR_2 PORT MAP (A => V1, B => V2, Y => Vv);
   VxorN : XOR_2 PORT MAP (A => Nn, B => Vv, Y => VxorredN);
	Zbarcase : AND_2 PORT MAP (A => Zzbar, B => VxorredN, Y => L1);
	lessthan : OR_2 PORT MAP (A => L1, B => Zz, Y => Ll);
	notlessthan : INVERTER PORT MAP (A => Ll, Y => nL);
	Greaterthan : OR_2 PORT MAP (A => nL, B => Zz, Y => G);
	Retrieve_Z : AND_2 PORT MAP (A => Zz, B => Zz, Y => Z);
	Retrieve_N : AND_2 PORT MAP (A => Nn, B => Nn, Y => N);
	Retrieve_V : AND_2 PORT MAP (A => Vv, B => Vv, Y => V);
	Retrieve_L : AND_2 PORT MAP (A => Ll, B => Ll, Y => L);
  end Struct;