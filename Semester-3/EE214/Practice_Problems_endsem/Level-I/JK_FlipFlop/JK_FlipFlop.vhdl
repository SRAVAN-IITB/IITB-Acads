library ieee;
use ieee.std_logic_1164.all;

entity JK_FlipFlop is
	port (J, K : in std_logic;
			Clock : in std_logic;
			Enable : in std_logic;
			Preset, Reset : in std_logic;
			Q, QBar : out std_logic);
end entity JK_FlipFlop;

architecture Behav of JK_FlipFlop is
signal temp : std_logic:= '0';
begin
	
	Q_CHANGE: process (clock, enable, preset, reset)
		
	begin
		if (Enable = '1') then
			if (clock =  '0' and clock'event) then
				if J = '1' then
					if K = '1' then
						temp <= not temp;
					else
						temp <= '1';
					end if;
				else
					if K='1' then
						temp <= '0';
					else
						null;
					end if;
				end if;
			end if;
--				case JK is
--					when "00" => null;
--					when "01" => temp <= '0';
--					when "10" => temp <= '1';
--					when "11" => temp <= not temp;
--				end case;
--				[This approach seems to cover only 4 out of 81 cases, which I am yet to figure out]
			if (preset = '1') then
					temp <= '1';
				elsif (reset = '1') then
					temp <= '0';
			end if;
		elsif (Enable = '0') then
			NULL;
		end if;

	end process Q_CHANGE;

	Q <= temp;
	QBAR <= not temp;
	
end Behav;
