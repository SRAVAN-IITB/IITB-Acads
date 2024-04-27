-----------------------------------------------------
-- Behvaior Design of GCD calculator
-- by Weijun Zhang, 05/2001
-----------------------------------------------------
library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_arith.all;
use work.all;
-----------------------------------------------------
entity gcd1 is
port(   Data_X:    in unsigned(3 downto 0);
        Data_Y:    in unsigned(3 downto 0);
        Data_out:  out unsigned(3 downto 0)
);
end gcd1;

architecture behv of gcd1 is
begin
    process(Data_X, Data_Y)
	variable tmp_X, tmp_Y: unsigned(3 downto 0);
    begin
	tmp_X := Data_X;
	tmp_Y := Data_Y;		 

	for i in 0 to 15 loop
	    if (tmp_X/=tmp_Y) then		
		if (tmp_X < tmp_Y) then
            	    tmp_Y := tmp_Y - tmp_X;
		else
		    tmp_X := tmp_X - tmp_Y;
		end if;
	    else
		Data_out <= tmp_X;
	    end if;
	end loop;
		
    end process;
end behv;
-----------------------------------------------------
