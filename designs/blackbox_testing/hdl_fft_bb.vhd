-- Generated from Simulink block hdl_fft_bb/Wideband_FFT_slim
library IEEE;
use IEEE.std_logic_1164.all;
library xil_defaultlib;
use xil_defaultlib.conv_pkg.all;
entity hdl_fft_bb_ip_struct is
  port (
    rst : in std_logic_vector( 1-1 downto 0 );
    in_sync : in std_logic_vector( 1-1 downto 0 );
    in_valid : in std_logic_vector( 1-1 downto 0 );
    shiftreg : in std_logic_vector(32-1 downto 0);
    in_im_0 : in std_logic_vector( 8-1 downto 0 );
    in_re_0 : in std_logic_vector( 8-1 downto 0 );
    clk_1 : in std_logic;
    ce_1 : in std_logic;
    out_sync : out std_logic;
    out_valid : out std_logic;
    ovflw : out std_logic;
    out_im_0 : out std_logic_vector( 16-1 downto 0 );
    out_re_0 : out std_logic_vector( 16-1 downto 0 )
  );
end hdl_fft_bb_ip_struct;
architecture structural of hdl_fft_bb_ip_struct is 
  component hdl_fft_bb_ip
    port (
        rst : in std_logic_vector( 1-1 downto 0 );
        in_sync : in std_logic_vector( 1-1 downto 0 );
        in_valid : in std_logic_vector( 1-1 downto 0 );
        shiftreg : in std_logic_vector(32-1 downto 0);
        in_im_0 : in std_logic_vector( 8-1 downto 0 );
        in_re_0 : in std_logic_vector( 8-1 downto 0 );
        clk_1 : in std_logic;
        ce_1 : in std_logic;
        out_sync : out std_logic;
        out_valid : out std_logic;
        ovflw : out std_logic; 
        out_im_0 : out std_logic_vector( 16-1 downto 0 );
        out_re_0 : out std_logic_vector( 16-1 downto 0 )
       ); 
    end component;
begin
  hdl_fft_bb_ip_inst : hdl_fft_bb_ip
  port map (
    rst => rst, 
    in_sync => in_sync, 
    in_valid => in_valid,
    shiftreg  => shiftreg,
    in_im_0 => in_im_0, 
    in_re_0 => in_re_0,
    clk => clk_l,
    out_sync => out_sync,
    out_valid => out_valid,
    ovflw => ovflw,
    out_im_0 => out_im_0,
    out_re_0 => out_re_0
  );
end structural;
