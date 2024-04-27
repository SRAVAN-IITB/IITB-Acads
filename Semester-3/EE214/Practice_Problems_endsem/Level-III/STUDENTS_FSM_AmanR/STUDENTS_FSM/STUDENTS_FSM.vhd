library ieee;
use ieee.std_logic_1164.all;

library work;
use work.Gates.all;

entity STUDENTS_FSM is
	port ( A,B,C,D,E,reset, clock : in std_logic; output : out std_logic);
end STUDENTS_FSM;

architecture STRUCT of STUDENTS_FSM is
	component DFF_JK is
		port( D, clock, enable, preset, reset : in std_logic; Q,QBAR : out std_logic);
	end component DFF_JK;

	component AND_9 is
		port (A,B,C,D,E,F,G,H,I : in std_logic; Y : out std_logic);
	end component AND_9;

	signal SBAR : std_logic_vector(3 downto 0);
	signal S3,S2,S1,S0,D0,D1,D2,D3 : std_logic;
	signal S3_CHANGE,S2_CHANGE, S1_CHANGE, S0_CHANGE : std_logic;
	signal AB,BB,CB,DB,EB,T1,T2,T3,T4,T5,T6,T7,T8,T9,O1,O2,O3 : std_logic;
	
begin
	
	D_3: DFF_JK port map ( D3, clock, '1', '0', '0',  S3, SBAR(3));
	D_2: DFF_JK port map ( D2, clock, '1', '0', '0',  S2, SBAR(2));
	D_1: DFF_JK port map ( D1, clock, '1', '0', '0',  S1, SBAR(1));
	D_0: DFF_JK port map ( D0, clock, '1', '0', '0',  S0, SBAR(0));

	AND9_1 :  AND_9 port map(SBAR(3),SBAR(2),	SBAR(1),	SBAR(0),		A,		BB, CB,	D, 	E, 	T1);
	AND9_2 :  AND_9 port map(SBAR(3),SBAR(2),	SBAR(1),	S0,			A,		BB, C, 	DB,	EB, 	T2);
	AND9_3 :  AND_9 port map(SBAR(3),SBAR(2),	S1,		S0,			A,		BB, C,	DB, 	E, 	T3);
	AND9_4 :  AND_9 port map(SBAR(3),SBAR(2),	S1,		SBAR(0),		AB,	BB, C, 	DB, 	EB, 	T4);
	AND9_5 :  AND_9 port map(SBAR(3),S2,		S1,		SBAR(0),		AB,	BB, C, 	DB, 	E, 	T5);
	AND9_6 :  AND_9 port map(SBAR(3),S2,		S1,		S0,			AB,	B,	 C, 	D, 	EB,	T6);
	AND9_7 :  AND_9 port map(SBAR(3),S2,		SBAR(1),	S0,			A,		BB, C,	DB, 	EB, 	T7);
	AND9_8 :  AND_9 port map(SBAR(3),S2,		SBAR(1),	SBAR(0),		A,		BB, CB,	D, 	E, 	T8);
	
	T9 <= (S3 and (SBAR(2) and (SBAR(1) and SBAR(0))));
	
	INV1: INVERTER port map(A,AB);
	INV2: INVERTER port map(B,BB);
	INV3: INVERTER port map(C,CB);
	INV4: INVERTER port map(D,DB);
	INV5: INVERTER port map(E,EB);
	
	OR0 : OR_2 port map(T9,T8,S3_CHANGE);
	OR1 : OR_2 port map(T4,T8,S2_CHANGE);
	OR2 : OR_2 port map(T6,T2,S1_CHANGE);
	OR3 : OR_2 port map(T1,T3,O1);
	OR4 : OR_2 port map(T5,T7,O2);
	OR5 : OR_2 port map(T9,O1,O3);
	OR6 : OR_2 port map(O2,O3,S0_CHANGE);
	
	XOR3: XOR_2 port map(S3_CHANGE,S3,D3);
	XOR2: XOR_2 port map(S2_CHANGE,S2,D2);
	XOR1: XOR_2 port map(S1_CHANGE,S1,D1);
	XOR0: XOR_2 port map(S0_CHANGE,S0,D0);
	
	OUTPUT <= S3;
	
end STRUCT;