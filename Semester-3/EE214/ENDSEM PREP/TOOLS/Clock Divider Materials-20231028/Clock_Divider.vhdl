library ieee;
use ieee.std_logic_1164.all;

entity Clock_Divider is
	port ( Clk_50  : in std_logic;
		   Reset   : in std_logic;
		   Clk_Out : out std_logic );
end Clock_Divider;

architecture Struct of Clock_Divider is
	
	begin
		CLKDIV: process(Clk_50)
					variable count : integer;
					variable outp  : std_logic;
					begin
						if (Clk_50'event and Clk_50 = '1') then  -- Clock is always the master input
																 --	All the operation is to be written inside Clk'event
							if (Reset = '1') then
								count := 0;			-- By reset it is advisable to put every value to 0 (Reset  State)
								outp  := '0';
							else
								if (count = 100000) then -- Change the count value according to given frequency (here f = 250Hz)
									outp  := not outp;
									count := 0;
								else
									outp  := outp;
									count := count + 1;
								end if;
							end if;
						end if;
						Clk_Out <= outp;
				end process;
end Struct;