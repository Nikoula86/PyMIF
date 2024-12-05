import glob
import xml.etree.ElementTree as et

def import_metadata_luxendo(
                 path: str,
                 ):
    """_summary_

    Args:
        path (str): _description_

    Returns:
        _type_: _description_
    """
    sample_meta = {}
    
    xml_file = glob.glob(f"{path}/*_bdv.xml")[0]
    xtree = et.parse(xml_file)
    xroot = xtree.getroot()
    # print(xroot)
    # print(xroot.tag)
    
    view_setups = xroot.find("SequenceDescription").find("ViewSetups")
    timepoints = xroot.find("SequenceDescription").find("Timepoints")
    voxel_info = view_setups.find("ViewSetup").find("voxelSize")
    
    sample_meta["size_t"] = int(timepoints.find("last").text)+1
    sample_meta["dt"] = 1.
    sample_meta["unit_t"] = "s"
    
    attributes = view_setups.findall("Attributes")
    attributes_names = [a.attrib["name"] for a in attributes]
    channel_attrib = attributes[attributes_names.index("channel")].findall("Channel")
    sample_meta["size_c"] = len(channel_attrib)
    for i, ch in enumerate(channel_attrib):
        sample_meta[f"channel_{i}_ID"] = ch.find("id").text
        sample_meta[f"channel_{i}_name"] = ch.find("name").text
        sample_meta[f"channel_{i}_color"] = str(i)
        
    xyz_size = view_setups.find("ViewSetup").find("size").text.split()
    d_xyz = voxel_info.find("size").text.split() 
    
    sample_meta["size_z"] = int(xyz_size[2])
    sample_meta["dz"] = float(d_xyz[2])
    sample_meta["unit_z"] = voxel_info.find("unit").text
    
    sample_meta["size_y"] = int(xyz_size[1])
    sample_meta["dy"] = float(d_xyz[1])
    sample_meta["unit_y"] = voxel_info.find("unit").text
    
    sample_meta["size_x"] = int(xyz_size[0])
    sample_meta["dx"] = float(d_xyz[0])
    sample_meta["unit_x"] = voxel_info.find("unit").text
    
    filename = sorted(glob.glob(f"{path}/*.lux.h5"))[0].split("/")[-1]
    filename = filename.replace("tp-0", "tp-{t}")
    filename = filename.replace("ch-0", "ch-{c}")
    sample_meta["file_template"] = filename

    return sample_meta
