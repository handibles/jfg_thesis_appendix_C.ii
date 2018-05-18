

cd /home/user/ad_comm/LEfSe/

mkdir ad_comm_test

cd ad_comm_test

python2.7 /home/user/ad_comm/LEfSe/format_input.py /user/ad_comm/ad_comm_lefse_input.txt FORMATTED.in -c 1 -s 2 -u 3 -o 1000000

python2.7 /home/user/ad_comm/LEfSe/run_lefse.py /home/user/ad_comm/LEfSe/ad_comm_test/FORMATTED.in /home/user/ad_comm/LEfSe/ad_comm_test/RESULT.res -l 3 --min_c 2

python2.7 /home/user/ad_comm/LEfSe/plot_res_mod.py /home/user/ad_comm/LEfSe/ad_comm_test/RESULT.res /home/user/ad_comm/LEfSe/ad_comm_test/IMAGE.png --dpi 300

python2.7 /home/user/ad_comm/LEfSe/plot_cladogram_mod.py /home/user/ad_comm/LEfSe/ad_comm_test/RESULT.res /home/user/ad_comm/LEfSe/ad_comm_test/RESULT.IMAGE.png --format png --dpi 300 --all_feats ALL_FEATS

