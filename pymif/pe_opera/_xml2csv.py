import pandas as pd
import xml.etree.ElementTree as et
import tqdm
import os


def xml2csv(exp_folder):

    xtree = et.parse(os.path.join(exp_folder, "Images", "Index.idx.xml"))
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

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ChannelID")
        row["channel"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}TimepointID")
        row["timepoint"] = int(x.text)

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ChannelName")
        row["chName"] = x.text

        x = image.find("{http://www.perkinelmer.com/PEHH/HarmonyV5}ExposureTime")
        row["expTime"] = float(x.text)

        df = pd.concat([df, pd.Series(row).to_frame().T], ignore_index=True)


    # print(df.head())
    df.to_csv(os.path.join(exp_folder, "metadata.csv"))

    # print("cioa")

#####################@###############################################################

if __name__ == '__main__':
    #####################

    ### mac gio
    # path = "/Volumes/sharpe/data/Vascular_micromass/Opera/TIMELAPSE/" "Timelapse4_041021/"
    # folder_raw = os.path.join(path)

    ### windows nicola
    path = os.path.join('data','Vascular_micromass','Opera','TIMELAPSE','Timelapse4_041021')
    folder_raw = os.path.join("X:", os.sep, path)

    exp_folder = os.path.join(
        "gio_Pecam-Sox9_20x-24h_041021__2021-10-04T16_06_44-Measurement_1"
    )

    # print(folder_raw)
    # print(exp_folder)

    #####################

