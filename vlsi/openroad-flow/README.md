# installing opendroad-flow 

# Installing Open PDK 

follow the [guide](http://opencircuitdesign.com/open_pdks/index.html)

# xschem sky130 integration 

follow the [guide](https://xschem.sourceforge.io/stefan/xschem_man/tutorial_xschem_sky130.html)

# Tutorial CMOS inverter using Magic and Sky130

[tutorial](http://web02.gonzaga.edu/faculty/talarico/vlsi/13_magic_inverter_sky130.pdf)

# installing OSS CAD Suite (yosys)

1. Download an archive matching your OS from [the releases page](https://github.com/YosysHQ/oss-cad-suite-build/releases/latest).
2. Extract the archive to a location of your choice (for Windows it is recommended that path does not contain spaces)
3. On macOS to allow execution of quarantined files ```xattr -d com.apple.quarantine oss-cad-suite-darwin-x64-yyymmdd.tgz``` on downloaded file, or run: ```./activate``` in extracted location once.
4. Set the environment as described below.

```
export PATH="<extracted_location>/oss-cad-suite/bin:$PATH"
```

# Working with OpenROAD-flow-scripts

```
cd OpenROAD-flow-scripts
source env.sh
export OPENROAD_EXE=$(command -v openroad)
export YOSYS_CMD=$(command -v yosys)
```
working example \w intro_to_vlsi/opendroad
```
cd OpenROAD-flow-scripts/flow
ESIGN_CONFIG=../../intro_to_vlsi/openroad/config.mk make
ESIGN_CONFIG=../../intro_to_vlsi/openroad/config.mk make gui_final
```

# Adi Teman lectures 
[youtube](https://www.youtube.com/watch?v=tywTQA_ko64&list=PLZU5hLL_713wQgIjRekOueTMJAyaoze3F&index=2)


