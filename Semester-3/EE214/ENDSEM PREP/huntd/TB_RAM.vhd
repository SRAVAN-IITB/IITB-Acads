--------------------------------------------------------------------
-- Test Bench for memory module (ESD book Chapter 5)
-- by Weijun Zhang, 04/2001
--
-- use loop statement to test module completely
-------------------------------------------------------------------- 

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.std_logic_arith.all;
use ieee.std_logic_unsigned.all;
 
entity MEM_TB is 		        -- entity declaration
end MEM_TB; 

--------------------------------------------------------------------
 
architecture TB of MEM_TB is 

component SRAM is	
port(	Clock:		in std_logic;	
	Enable:		in std_logic;
	Read:		in std_logic;
	Write:		in std_logic;
	Read_Addr:	in std_logic_vector(1 downto 0);
	Write_Addr: 	in std_logic_vector(1 downto 0);			 
	Data_in: 	in std_logic_vector(3 downto 0);
	Data_out: 	out std_logic_vector(3 downto 0)
);	
end component;

signal T_Clock, T_Enable, T_Read, T_Write: std_logic;
signal T_Data_in, T_Data_out: std_logic_vector(3 downto 0);
signal T_Read_Addr: std_logic_vector(1 downto 0);
signal T_Write_Addr: std_logic_vector(1 downto 0);

begin 
	
    U_CKT: SRAM port map (T_Clock, T_Enable, T_Read, T_Write,
		T_Read_Addr, T_Write_Addr, T_Data_in, T_Data_out);

    Clk_sig: process
    begin
        T_Clock<='1';		        -- clock cycle 10 ns
        wait for 5 ns;
        T_Clock<='0';
        wait for 5 ns;
    end process;
							
    process 
        variable err_cnt: integer := 0;
    begin
    	
	T_Enable <= '1';   
	T_Read <= '0';
	T_Write <= '0';
	T_Write_Addr <= (T_Write_Addr'range => '0');
	T_Read_Addr <= (T_Read_Addr'range => '0');
	T_Data_in <= (T_Data_in'range => '0');		
	wait for 20 ns;
	
	-- test write		
	for i in 0 to 3 loop
	    T_Write_Addr <= T_Write_Addr + '1';
	    T_Data_in <= T_Data_in + "10";
	    T_Write <= '1';
	    wait for 10 ns;				
	    assert (T_Data_out="ZZZZ") 
                report "Something wrong!" severity Error;
            if (T_Data_out /= "ZZZZ") then
                err_cnt := err_cnt + 1;
            end if; 
	end loop;
		
	-- test read
	for i in 0 to 2 loop
	    T_Read_Addr <= T_Read_Addr + '1';
	    T_Read <= '1';
	    wait for 10 ns;	
	    assert (conv_integer(T_Data_out)=2*conv_integer(T_Read_Addr))
		report "Something wrong!" severity Error;
	    if (conv_integer(T_Data_out)/=2*conv_integer(T_Read_Addr)) then
		err_cnt := err_cnt + 1;
	    end if;			  
	end loop;
		
	-- summary of all the tests
        if (err_cnt=0) then                     
            assert false 
            report "Testbench of ROM completed successfully!" 
            severity note; 
        else
            assert true
            report "Something wrong, try again"
            severity error;
        end if;

	wait;

    end process;

end TB;

--------------------------------------------------------------------------
configuration CFG_TB of MEM_TB is 
        for TB 
        end for; 
end CFG_TB; 
--------------------------------------------------------------------------



