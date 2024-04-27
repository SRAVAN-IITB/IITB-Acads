library IEEE;
use IEEE.std_logic_1164.all;

entity JK_FF is
	port(clock, J, K, Reset, Preset, Enable : in std_logic; Q, QBar : out std_logic);
end JK_FF;

architecture BHV of JK_FF is
	signal buff, outp : std_logic;
begin
	Data_change : process(clock, Enable)
	begin
		
--		if(Preset <= '1') then
--		buff <= '1';
--		end if;
--		
--		if(Reset <= '1') then
--		buff <= '0';
--		end if;
		
		if(Enable = '1') then
			if(clock = '0' and clock'event) then
				if(J = '0' and K = '0') then
					NULL;
				end if;
				
				if(J = '0' and K = '1') then
					buff <= '0';
				end if;
				
				if(J = '1' and K = '0') then
					buff <= '1';
				end if;
				
				if(J = '1' and K = '1') then
					buff <= not(buff);
				end if;
				
			end if;
		end if;
	
	end process Data_change;
		
	outp <= (buff and (not Reset)) or Preset;
		Q <= outp;
	QBar <= not(outp);
	
end BHV;