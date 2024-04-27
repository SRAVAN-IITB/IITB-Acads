library IEEE;
use IEEE.std_logic_1164.all;

entity Testbench is
end entity;

architecture TB_BHV of Testbench is
component D_FlipFlop is
	port (D : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Reset, Preset : in std_logic;
			Q, QBar : out std_logic);
end component D_FlipFlop;

signal Reset, Preset, D, Enable, Clock : std_logic := '0';
signal Output : std_logic_vector(1 downto 0);

begin	
	DFF_Instance : D_FlipFlop port map (D => D, 
													Clock => Clock,
													Reset => Reset, Preset => Preset,
													Enable => Enable,
													Q => Output(0), QBar => Output(1));
	
		clock <= not clock after 10ns;
		D <= not D after 8ns;
		Preset <= '1', '0' after 15ns;
		Reset <= '1', '0' after 30ns;
		Enable <= '1';
	
	end_sim : process
	begin
		wait for 400ns;
		ASSERT False
		REPORT "Sim ended"
		SEVERITY failure;
	end process end_sim;

end TB_BHV;