library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end Testbench;

architecture TB_BHV of Testbench is
	component Synch_UpDown_Ctr_16bit is
		PORT(Clock, Enable, Preset, Reset, M : in std_logic;
			Counter : out std_logic_vector(15 downto 0));
	end component Synch_UpDown_Ctr_16bit;
	
	signal Clock, Preset , Reset, M : std_logic := '0';
	signal Enable : std_logic;
	signal Counter : std_logic_vector(15 downto 0);
	
	begin
		counterInst : Synch_UpDown_Ctr_16bit PORT MAP(Clock => Clock, Enable => Enable,
																	 Preset => Preset, Reset => Reset,
																	 M => M, Counter => Counter);
		
		clock  <= not clock after 10ns;
		Enable <= '1';
		M <= '0', '1' after 1310720ns;
	
--		Counter <= "0000000000000000", "1111111111111111" after 400ns;
		
		end_sim : process
			 begin
				wait for 2621440ns;
				ASSERT false
				REPORT "Sim ended"
				SEVERITY Failure;
		end process end_sim;
	
end TB_BHV;