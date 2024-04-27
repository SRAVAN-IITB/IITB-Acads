library IEEE;
use IEEE.std_logic_1164.all;
use work.gates.all;

entity Synch_Down_Ctr_4bit is
	PORT(clock : in std_logic;
			counter : out std_logic_vector(3 DOWNTO 0));
end Synch_Down_Ctr_4bit;

architecture STRUCT of Synch_Down_Ctr_4bit is
	component D_FF is
		PORT(clock, D, DBar, Reset, Preset, Enable : in std_logic; Q, QBar : out std_logic);
	end component D_FF;

	signal Q 	: std_logic_vector(3 DOWNTO 0) := "1111";
	signal QBar : std_logic_vector(3 DOWNTO 0) := "0000";
	signal M : std_logic := '0';									-- M = 0 ensures down-counting

	begin
		D_FF1 : D_FF PORT MAP(clock => clock, 
										D => '1', DBar => '1', 
										Reset => '0', Preset => '0', 
										Enable => '1', 
										Q => Q(0), QBar => QBar(0));
										
--		D_FF2 : D_FF PORT MAP(clock => clock, 
--										D => (Q(0) AND M) or (QBar(0) AND (NOT(M))), 
--										DBar => (Q(0) AND M) or (QBar(0) AND (NOT(M))), 
--										Reset => '0', Preset => '0', 
--										Enable => '1', 
--										Q => Q(1), QBar => QBar(1));
--										
--		D_FF3 : D_FF PORT MAP(clock => clock, 
--										D => (Q(1) AND Q(0) AND M) or (QBar(1) AND QBar(0) AND (NOT(M))), 
--										DBar => (Q(1) AND Q(0) AND M) or (QBar(1) AND QBar(0) AND (NOT(M))),
--										Reset => '0', Preset => '0', 
--										Enable => '1', 
--										Q => Q(2), QBar => QBar(2));
--										
--		D_FF4 : D_FF PORT MAP(clock => clock, 
--										D => (Q(2) AND Q(1) AND Q(0) AND M) or (QBar(2) AND QBar(1) AND QBar(0) AND (NOT(M))), 
--										DBar => (Q(2) AND Q(1) AND Q(0) AND M) or (QBar(2) AND QBar(1) AND QBar(0) AND (NOT(M))), 
--										Reset => '0', Preset => '0', 
--										Enable => '1', 
--										Q => Q(3), QBar => QBar(3));
		
		counter <= Q;

end STRUCT;
