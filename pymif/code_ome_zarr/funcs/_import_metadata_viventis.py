import glob
import xml.etree.ElementTree as et

def import_metadata_viventis(
                            path: str,
                            ):
    """_summary_

    Args:
        path (str): _description_

    Returns:
        _type_: _description_
    """
    sample_meta = {}
    
    xml_file = glob.glob(f"{path}/*.companion.ome")[0]
    xtree = et.parse(xml_file)
    xroot = xtree.getroot()
    # print(xroot.tag)
    # print(xroot.attrib)
    
    for child in xroot:
        for cchild in child:
            # print(cchild.tag, cchild.attrib)
            d = cchild.attrib
            sample_meta["size_t"] = int(d["SizeT"])
            sample_meta["size_c"] = int(d["SizeC"])
            sample_meta["size_z"] = int(d["SizeZ"])
            sample_meta["size_y"] = int(d["SizeY"])
            sample_meta["size_x"] = int(d["SizeX"])
            sample_meta["dt"] = float(d["TimeIncrement"])
            sample_meta["dz"] = float(d["PhysicalSizeZ"])
            sample_meta["dy"] = float(d["PhysicalSizeY"])
            sample_meta["dx"] = float(d["PhysicalSizeX"])
            sample_meta["unit_t"] = d["TimeIncrementUnit"]
            sample_meta["unit_z"] = d["PhysicalSizeZUnit"].replace("µ","u")
            sample_meta["unit_y"] = d["PhysicalSizeYUnit"].replace("µ","u")
            sample_meta["unit_x"] = d["PhysicalSizeXUnit"].replace("µ","u")
            i=0
            for ccchild in cchild:
                if "Channel" in ccchild.tag:
                    d = ccchild.attrib
                    sample_meta[f"channel_{i}_ID"] = d["ID"]
                    sample_meta[f"channel_{i}_name"] = d["Name"]
                    sample_meta[f"channel_{i}_color"] = str(d["Color"])
                    # print(d.tag, d.attrib)
                    i+=1
                    
    xmlns = xroot.tag[:-3]
    filename = xroot.find(f"{xmlns}Image").find(f"{xmlns}Pixels").find(f"{xmlns}TiffData").find(f"{xmlns}UUID").attrib["FileName"]
    filename = filename.replace("t0001", "t{t}")
    filename = filename.replace(sample_meta["channel_0_name"], "{c}")
    # print(xroot.find("Image").find("Pixel").find("TiffData").find("UUID").attrib)
    sample_meta["file_template"] = filename
        
    return sample_meta