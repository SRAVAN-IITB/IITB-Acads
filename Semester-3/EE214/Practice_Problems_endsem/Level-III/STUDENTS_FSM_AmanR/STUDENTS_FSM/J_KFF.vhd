library ieee;
use ieee.std_logic_1164.all;

entity J_KFF is 
	port(clock, J , K, enable, preset, reset : in std_logic; Q , QBAR : out std_logic);
end J_KFF;

architecture BHV of J_KFF is

signal temp : std_logic:= '0';
begin
	
	Q_CHANGE: process (clock,enable,preset,reset)
		
	begin
		if (enable = '1') then
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
		end if;
		if (preset = '1') then
				temp <= '1';
			elsif (reset = '1') then
				temp <= '0';
		end if;
	end process Q_CHANGE;

	Q <= temp;
	QBAR <= not temp;
	
end BHV;