library IEEE;
use IEEE.std_logic_1164.all;
library work;
use work.gates.all;

entity Synch_UpDown_Ctr_16bit is
	PORT(Clock, Enable, Preset, Reset, M : in std_logic;
			Counter : out std_logic_vector(15 downto 0));
end entity;

architecture STRUCT of Synch_UpDown_Ctr_16bit is
	component JK_FlipFlop is
		port (J, K : in std_logic;
				Clock : in std_logic;
				Enable : in std_logic;
				Preset, Reset : in std_logic;
				Q, QBar : out std_logic);
	end component JK_FlipFlop;

	signal Mbar_OR_Q, M_OR_QBar, upORdown : std_logic_vector(15 downto 0);
	signal    Qout : std_logic_vector(15 downto 0) := "0000000000000000";
	signal QBarout : std_logic_vector(15 downto 0) := "1111111111111111";
	--	The trouble is that I don't know how to initialise counter up-count and down-count. 
	--	As of now I'm trying for atleast one way to work but simulation shows undefined output :(
	
	begin
		JK0 : JK_FlipFlop PORT MAP(J => Enable, K => Enable,
											Clock => Clock,
											Enable => Enable,
											Preset => Preset, 
											Reset => Reset,
											Q => Qout(0), QBar => QBarout(0));
		
		Mbar_OR_Q(0) <= not(M) and Qout(0);
		M_OR_QBar(0) <= M and QBarout(0);
		upORdown(0)	 <= (not(M) and Qout(0)) or (M and QBarout(0));
		
		JK : for i in 1 to 15 generate
					andforQ : AND_2 PORT MAP(A => Mbar_OR_Q(i-1), B => Qout(i-1), 
													 Y => Mbar_OR_Q(i));
													 
					andforQBar : AND_2 PORT MAP(A => M_OR_QBar(i-1), B => QBarout(i-1), 
													 Y => M_OR_QBar(i));
					
					ORfeedingJK : OR_2 PORT MAP(A => Mbar_OR_Q(i), B => M_OR_QBar(i),
													 Y => upORdown(i));
													 
					JK : JK_FlipFlop PORT MAP(J => upORdown(i), 
													  K => upORdown(i),
														Clock => Clock,
														Enable => Enable,
														Preset => '0', 
														Reset => '0',
														Q => Qout(i), 
														QBar => QBarout(i));
		end generate;
		
		Counter <= Qout;
end architecture STRUCT;