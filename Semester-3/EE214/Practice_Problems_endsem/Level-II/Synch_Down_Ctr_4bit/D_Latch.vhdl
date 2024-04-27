library ieee;
use ieee.std_logic_1164.all;

entity D_Latch is 
	port (clock, D :in std_logic ; Q, QBAR : out std_logic);
end D_Latch;

architecture BHV of D_Latch is
	signal buff : std_logic;
begin
	Data_Change: process (clock)
	begin
--		b <= c;
--		a <= b;
--		b <= d;
		if (clock = '0' and clock'event) then
			buff <= D;
		end if;	
	end process Data_Change;

	Q <= buff;
	QBAR <= not buff;

end bhv;