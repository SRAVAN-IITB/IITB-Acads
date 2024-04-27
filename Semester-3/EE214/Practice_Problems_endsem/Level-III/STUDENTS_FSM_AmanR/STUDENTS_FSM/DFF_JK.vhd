library ieee;
use ieee.std_logic_1164.all;

library work;
use work.Gates.all;

entity DFF_JK is
	port( D, clock, enable, preset, reset : in std_logic; Q,QBAR : out std_logic);
end DFF_JK;

architecture BHV of DFF_JK is
	component J_KFF is 
		port(clock, J , K, enable, preset, reset : in std_logic; Q , QBAR : out std_logic);
	end component J_KFF;
	signal R,RBAR,RESETBAR, K_IN, O1, TEMP : std_logic;
begin
	
	INV1 : INVERTER port map(D, K_IN);
	JKFF_2: J_KFF port map (clock, D, K_IN, enable, '0', '0', R, RBAR);
	INV3 : INVERTER port map(reset, RESETBAR);

	AND1: AND_2 port map(RESETBAR, R, O1);
	OR1: OR_2 port map(preset, O1, TEMP);
	INV2 : INVERTER port map(TEMP, QBAR);
	Q <= TEMP;
end BHV;