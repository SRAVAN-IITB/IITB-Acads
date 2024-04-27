library IEEE;
use IEEE.std_logic_1164.all;

entity Asynch_Up_Ctr_4bit is
	port (Clock, Reset : in std_logic;
			Counter : out std_logic_vector(3 downto 0));
end Asynch_Up_Ctr_4bit;

architecture STRUCT of Asynch_Up_Ctr_4bit is
	component JK_FF is
		port(clock, J, K, Reset, Preset, Enable : in std_logic; Q, QBar : out std_logic);
	end component JK_FF;
--	
	signal Q_out : std_logic_vector(3 downto 0);
	signal QBar_out : std_logic_vector(3 downto 0);	
	
	begin	
	JKFF1 : JK_FF port map (J => '1', K => '1',
											Clock => Clock,
											Enable => '1',
											Preset => '0',
											Reset => Reset,
											Q => Q_out(0), QBar => QBar_out(0));
											
--	GenFF : for i in 1 to 3 generate
--	JKFF_Inst : JK_FF port map (J => '1', K => '1',
--												 Clock => Q_out(i-1),
--												 Enable => '1',
--												 Preset => '0',
--												 Reset => Reset,
--												 Q => Q_out(i), QBar => QBar_out(i));
--	end generate;
--	
	
	Counter <= Q;
	
end STRUCT;