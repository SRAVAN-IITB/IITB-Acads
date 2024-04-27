---------------------------------------------------
-- 3-bit Shift-Register/Shifter
-- (ESD book figure 2.6)
-- by Weijun Zhang, 04/2001
--
-- reset is ignored according to the figure
---------------------------------------------------

library ieee ;
use ieee.std_logic_1164.all;

---------------------------------------------------

entity shift_reg is
port(	I:		in std_logic;
	clock:		in std_logic;
	shift:		in std_logic;
	Q:		out std_logic
);
end shift_reg;

---------------------------------------------------

architecture behv of shift_reg is

    -- initialize the declared signal
    signal S: std_logic_vector(2 downto 0):="111";

begin
    
    process(I, clock, shift, S)
    begin

	-- everything happens upon the clock changing
	if clock'event and clock='1' then
	    if shift = '1' then
		S <= I & S(2 downto 1);
	    end if;
	end if;

    end process;
	
    -- concurrent assignment
    Q <= S(0);

end behv;

----------------------------------------------------
