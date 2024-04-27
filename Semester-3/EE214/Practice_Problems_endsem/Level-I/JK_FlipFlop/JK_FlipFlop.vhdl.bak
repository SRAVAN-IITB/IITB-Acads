library IEEE;
use IEEE.std_logic_1164.all;

entity JK_FlipFlop is
	port (J, K : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Preset, Reset : in std_logic;
			Q, NotQ : out std_logic);
end entity JK_FlipFlop;

architecture BHV of JK_FlipFlop is
signal NextState, PrevState : std_logic;
signal Input : std_logic_vector(1 downto 0);
begin
	Input <= J & K;
	
	FlipFlopProcess : process (Clock, Enable, Preset, Reset, NextState)
	begin
		if (Preset = '1') then
			NextState <= '0';
		elsif (Reset = '1') then
			NextState <= '1';
			
		elsif (Clock = '1' and Clock'event) then
			if (Enable  = '1') then
				case (Input) is
					-- Set and Reset the FlipFlop
					when "10" =>
						NextState <= '1';
					when "01" => 
						NextState <= '0';
						
					-- Memory Storage Case
					when "00" => 
						NextState <= PrevState;
						
					-- Toggle Case
					when "11" =>
						NextState <= not PrevState;
					
					when others =>
						null;
				end case;
			end if;
		end if;
		
		Q <= NextState;
		NotQ <= not NextState;
		PrevState <= NextState;
	end process;
end BHV;
