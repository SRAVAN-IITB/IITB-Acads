------------------------------------------------------------------
-- Test Bench for 2-bit register (ESD figure 2.6)
-- by Weijun Zhang, 04/2001
--
-- several ways you may use to specify the concurrent clock signal
------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;
use ieee.std_logic_unsigned.all;
use ieee.std_logic_arith.all;

entity reg_TB is			-- entity declaration
end reg_TB;

------------------------------------------------------------------

architecture TB of reg_TB is

    component reg
    port(   I:		in std_logic_vector(1 downto 0);
    	    clock:	in std_logic;
	    load:	in std_logic;
	    clear:	in std_logic;
	    Q:		out std_logic_vector(1 downto 0)
    );
    end component;

    signal T_I:		std_logic_vector(1 downto 0);
    signal T_clock:	std_logic;
    signal T_load:	std_logic;
    signal T_clear:	std_logic;
    signal T_Q:		std_logic_vector(1 downto 0);
	
begin

    U_reg: reg port map (T_I, T_clock, T_load, T_clear, T_Q);
	
    -- concurrent process to offer the clock signal
    process
    begin
 	T_clock <= '0';
 	wait for 5 ns;
	T_clock <= '1';
	wait for 5 ns;
    end process;
	
    process							
							 
	variable err_cnt: integer :=0;
	
    begin								
	
	T_I <= "10";
	T_load <= '0';
	T_clear <= '1';
		
	-- case 1
	wait for 20 ns;
	T_load <= '1';
	wait for 10 ns;
	assert (T_Q="10") report "Test1 Failed!" severity error;
	if (T_Q/=T_I) then
	    err_cnt := err_cnt+1;
	end if;
		
	-- case 2				
	wait for 10 ns;
	T_load <= '0';
	wait for 10 ns;
	assert (T_Q="10") report "Test2 Failed!" severity error;
	if (T_Q/=T_I) then
	    err_cnt := err_cnt+1;
	end if;		
		
	-- case 3
	wait for 10 ns;
	T_clear <= '0';										   
	wait for 10 ns;
	assert (T_Q="00") report "Test3 Failed!" severity error;
	if (T_Q/=T_I) then
	    err_cnt := err_cnt+1;
	end if;
		
	-- case 4
	wait for 10 ns;
	T_clear <= '1';
	wait for 10 ns;
	assert (T_Q="00") report "Test4 Failed!" severity error;
	if (T_Q/=T_I) then
	    err_cnt := err_cnt+1;
	end if;
		
	-- summary of all the tests
	if (err_cnt=0) then	
	    assert false
	    report "Testbench of register completely successfully!"
	    severity note;
	else
	    assert true
	    report "Something wrong, check again pls!"
	    severity error;
	end if;
	
        wait;
		
    end process;

end TB;

------------------------------------------------------------------
configuration CFG_TB of reg_TB is
    for TB
    end for;
end CFG_TB;
------------------------------------------------------------------
