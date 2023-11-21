import pandas as pd
import xml.etree.ElementTree as et
import tqdm
import os, glob


def xml2csv(exp_folder,
            image_folder = "Images",
            meta_file_name = "metadata.csv",
            save = True):

    # print(os.path.join(exp_folder, image_folder, "*.xml"))
    xml_file = glob.glob(os.path.join(exp_folder, image_folder, "*.xml"))[0]
    xtree = et.parse(xml_file)
    xroot = xtree.getroot()

    images = xroot.findall("{http://www.perkinelmer.com/PEHH/HarmonyV5}Images")[0]
    print("Found %d images."%len(images))

    df = pd.DataFrame(
        {
            "filename": [],
            "Xpos": [],
            "Ypos": [],
            "Zpos": [],
            "row": [],
            "col": [],
            "field": [],
            "plane": [],
            "channel": [],
            "chName": [],
            "expTime": [],
        }
    )


    for i, image in tqdm.tqdm(enumerate(images.iter("{http://www.perkinelmer.com/PEHH/HarmonyV5}Image")), total=len(images)):
        # print(image.tag, image.attrib)

        row = {}
        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}URL")
        row["filename"] = x.text

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}PositionX")
        row["Xpos"] = float(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}PositionY")
        row["Ypos"] = float(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}PositionZ")
        row["Zpos"] = float(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}Row")
        row["row"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}Col")
        row["col"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}FieldID")
        row["field"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}PlaneID")
        row["plane"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}TimepointID")
        row["timepoint"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ChannelID")
        row["channel"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ChannelName")
        row["chName"] = x.text

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ChannelType")
        row["chType"] = x.text

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}MainExcitationWavelength")
        row["chWavelength"] = x.text

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ExposureTime")
        row["expTime"] = float(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ImageResolutionX")
        row["pixelSize"] = float(x.text)*1e6

        df = pd.concat([df, pd.Series(row).to_frame().T], ignore_index=True)


    # print(df.head())
    if save:
        df.to_csv(os.path.join(exp_folder, meta_file_name))

    return df

#####################################################################################

if __name__ == '__main__':
    #####################

    ### windows nicola
    # exp_folder = "/g/trivedi/Kristina_Stapornwongkul/ImageAnalysis/gastr_hcr_volumes/data/primary/date-20220304_hpa-96_plate-1_exp-1"
    exp_folder = "PATH-TO-EXPERIMENT"

    xml2csv(exp_folder,
        image_folder = "Images",
        meta_file_name = "metadata_PE.csv",
        save = True) 

    #####################

