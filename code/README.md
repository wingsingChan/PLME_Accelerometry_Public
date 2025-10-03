# Supplementary Information: Accelerometer and Data Logger for Small- to Medium-Sized Animal Research

## Overview
This document describes a low-cost, low-power, and customizable data logger designed for animal research. It primarily focuses on collecting movement data from small- to medium-sized animals using a compact accelerometer.

## Key Features
- **Light Weight**: 
  - Without battery or case: < 3 grams
  - With 290 mAh battery and case: 11 - 16 grams
  - With 500 mAh battery and case: 20 - 21 grams
- **Small Dimensions**: 
  - Without battery or case: 20mm x 20mm x 8mm
- **Extended Battery Performance**: Approximately 2 weeks to 3 weeks
- **Budget-friendly**: Less than $60 USD
- **User-friendly**: Made using easily purchased commercially-available parts from TinyCircuits, and can be customized using the Arduino ecosystem.

## Motivation
Understanding animal movements can provide valuable insights for biologists. Existing commercial accelerometers are often too expensive for widespread use, while DIY solutions may be bulky or have inadequate battery life for meaningful research. This project aims to overcome these limitations.

## Hardware Components
The accelerometer comprises four parts. All parts can be purchased from TinyCircuits or any of their retailers, except for the standard micro-SD card. The parts are:

1. **TinyCircuits TinyZero Processor Board with Accelerometer (ASM2021-R)**: An Arduino-compatible microcontroller with a SAMD21 processor and an integrated accelerometer. Note that there are two versions of this board available. It is important that you purchase the version of this board that includes the accelerometer, as this is what you will use to collect the data. 

2. **TinyCircuits MicroSD TinyShield (ASD2201-R)**: An adapter that allows you to conveniently save data to a MircoSD card.

3. **TinyCircuits Lithium Ion Polymer Battery**: TinyCircuits sells batteries with the appropriate connector on their website. You can purchase any size that fits your purpose, but there is a trade-off between battery capacity and weight. We found a good balance to be the 290mAh battery from TinyCircuits. It provided more than 2 weeks of battery life.

4. **Generic Micro SD Card**: For data storage.A 64GB card is advisable. 

#### Assembly
Once you have the parts, assembly is straight-forward. The MicroSD shield fits on top of the TinyZero processor board through the white connectors on each board. The battery plugs into the corresponding battery connector on the TinyZero processor board. The SD Card fits into the SD slot on the MicroSD shield.

#### Waterproofing
To make the logger suitable for field use, particularly on aquatic turtles, the assembly is waterproofed using:

- a small zip-lock bag,
- a thin layer of Parafilm, and 
- encased with using epoxy putty glue

## Software Setup
To operate the data logger, follow the below steps:

### 1. Prepare the Arduino Environment
- Install the Arduino Software IDE. 
- Open Arduino IDE. Go to `File` -> `Preferences`. 
- Copy and paste the following URL to the field `Additional Boards Manager URLs`, https://files.tinycircuits.com/ArduinoBoards/package_tinycircuits_index.json. 
- Click `OK`.
- Go to `Tools` -> `Board` -> `Boards Manager`. 
- Select and install **Arduino SAMD Boards (32-bits ARM Cortex-M0+)**. 
- On the same page, select and install **TinyCircuits SAMD Boards**. 

### 2. Add TinyZero Board Profiles into the Arduino IDE
- Connect the TinyZero Board to your computer using microUSB cable. Make sure you turn on the power switch on the side of your TinyZero. 
- Go to `Tools` -> `Board` -> `TinyZero`.
- Go to `Tools` -> `Port` and select the port labelled `COMXX (TinyScreen+)`. If it does not show up, ensure that your TinyZero is powered on, and you are using a microUSB cable that transmits data.

### 3. Install the Required Libraries
In the Aruino IDE, go to `Sketch` -> `Include Library` -> `Include .ZIP Library` and select the ZIP files you downloaded. Download and install each of the following libraries:

- **BMA250 Library**
- **Adafruit Sleepy Dog Library**
- **Arduino Low Power Library**
- **RTC Zero Library**

Make sure you installed the specific version of the library provided. You can find the ZIP files for these libraries in this project repository or on [Hackster.io](https://www.hackster.io/brian-k2/accelerometer-and-data-logger-for-small-animal-research-c877c6). Using other versions may lead to compatibility issues, which could result in the failure of the running model.

### 4. Upload the Software
- Connect the TinyZero board to a computer via a microUSB cable and select the correct board and port in the Arduino IDE (as described in step 2).
- Open **sdlogger.ino** in the Arduino editor. You can download the file in this project repository or on [Hackster.io](https://www.hackster.io/brian-k2/accelerometer-and-data-logger-for-small-animal-research-c877c6). 
- Click the "right arrow" symbol in the top menu to upload the code to the TinyZero board.

If the code does not upload, ensure that you have the correct board and port selected under `Tools`. If the issue persists, you may need to place the board into "Bootloader" mode. To do this, turn the power switch to the off position. Locate the small button on the bottom of the TinyZero board, press and hold the button as you move the power switch to the on position. The board should now be available in the Arduino editor.

### 5. Operating the Data Logger
Once the code is uploaded, you are ready to collect data! Whenever the power switch is in the "on" position, the accelerometer will either collect data or enter sleep mode based on the specified parameters.

To retrieve the data, turn the power switch off, then remove the microSD card and insert it into a computer with an microSD card reader. You should find a file named `logfile.csv` (or the name you designated in `filename` in the **sdlogger.ino**). This can be read with any program that supports CSV formats.

Note that if the data logger has operated for a duration shorter than the `SAVE_INTERVAL` parameter, the data will not be saved on the card.

### Configuring the Data Logger
The code provided in **sdlogger.ino** comprises several key sections:

- **Libraries**: Import necessary libraries for functioning.
- **Setup**: Initializes the SD card, RTC, and accelerometer.
- **Loop**: Collects data, writes it to the SD card, and manages sleep cycles.

#### Libraries
The first part of the code accesses the necessary libraries:

```cpp
#include <Wire.h>
#include <BMA250.h>
#include <RTCZero.h>            // Can only go into sleep on second intervals
                                // Does not have power down flash or SYSTICK modifications
#include <Adafruit_SleepyDog.h> // Has power down flash correct, but not SYSTICK
                                // Allows you to do millisecond sleep, but limited to WDT limit ~18s
#include <SdFat.h>
#include <SPI.h>
#include <ArduinoLowPower.h>    // Use RTC library for timing, therefore can only sleep in second intervals
                                // Has power down flash and SYSTICK correct
```

- **Wire.h** and **BMA250.h**: Used to access the accelerometer.
- **RTCZero.h**: Used to keep time with the Real Time Clock and to manage the initial sleep period.
- **Adafruit_SleepyDog.h**: Used to sleep between each accelerometer sample.
- **SdFat.h** and **SPI.h**: Used to access the SD card.
- **ArduinoLowPower.h**: Used to sleep between recording intervals.

You may notice that the RTCZero, Adafruit_SleepyDog, and ArduinoLowPower libraries all have sleep functions, but they differ in their implementation and usage. To make this code easy to use and accessible, we chose to stick with common libraries rather than create our own, utilizing each library where appropriate.

The **RTCZero library** uses the Real Time Clock and can only perform sleep periods in one-second intervals. Additionally, it has not fixed a well-known bug (see "SAMD21 systick" for more information) that can prevent the processor from waking up after sleeping. Therefore, we do not use it for sleeping, but rather to keep time with the Real Time Clock.

The **Adafruit_SleepDog library** utilizes the watchdog timer, allowing for much shorter sleep durations, which we use to sleep between each data point. However, it cannot sleep for longer than 18 seconds due to its reliance on the watchdog timer. It suffers from the same bug that has not been fixed in the source, so we include a modified version.

The **ArduinoLowPower library** does not have the well-known bugs present in the other libraries. However, like the RTCZero library, it uses the Real Time Clock for sleeping and can only perform sleep periods in one-second intervals. Thus, we only use it to sleep between recording intervals.

#### Parameters

```cpp
#define INITIAL_SLEEP_TIME 240 // hours; How long to sleep on power, allows to delay start of recording to save power
#define RECORD_INTERVAL 15 // seconds; How long to collect data for during each recording interval
#define SLEEP_BETWEEN_INTERVAL 900 // seconds; How long to sleep between recording intervals
#define SAVE_INTERVAL 14400 // seconds; How long to delay between SD card saves. Each save uses a lot of power, so we only save occasionally
#define SLEEP_BETWEEN_SAMPLES 60 // milliseconds; Time between samples, determined by data rate of accelerometer
#define TIME_BETWEEN_SAMPLES 64 // milliseconds; Adjustment factor to keep timestamp accurate after sleep
#define DATA_COUNT 4 // Number of data fields you are collecting. For x,y,z accelerations + temperature -> DATA_COUNT=4
```

To understand the parameters, it is helpful to know how the data logger works. It alternates between three activities: sleep, data collection, and saving data. The logger collects data in intervals with a length determined by the parameter `RECORD_INTERVAL`. After collecting data, it goes to sleep for a duration specified by `SLEEP_BETWEEN_INTERVAL`. It then wakes up and collects data again. To save power, it only saves data occasionally, after a duration equal to the parameter `SAVE_INTERVAL`.

Additionally, the first sleep period after the power switch is turned on can differ from all subsequent sleep periods. This is determined by the `INITIAL_SLEEP_TIME` parameter.

In the example above, once the power switch is turned on, the data logger sleeps for 240 hours (or 10 days) as specified by `INITIAL_SLEEP_TIME`. It then wakes up and collects data for 15 seconds, as determined by `RECORD_INTERVAL`. Afterward, it sleeps again for 900 seconds, as specified by `SLEEP_BETWEEN_INTERVAL`. This cycle of record-sleep-record-sleep continues as long as the battery remains functional. Every 4 hours (14400 seconds), it saves all the collected data to the SD card.

To change the recording and sleep behavior of the data logger, you can adjust the parameters above. Just be aware that more frequent recordings, longer recording intervals, and more frequent saves will consume more battery life.

#### Setup

```cpp
void setup() {
  // If SD initialization fails, light LED
  if (!sd.begin(CHIP_SELECT, SD_SCK_MHZ(50))) {
    SPI.end();
    lightLED();
  }
  // If file creation fails, light LED
  if (!file.open(fileName, O_WRONLY | O_CREAT | O_APPEND)) {
    SPI.end();
    lightLED();
  }

  detectReset(); // Checks reason for reset, for debugging
  file.sync();

  // Initialize RTC
  rtc.begin();
  rtc.setEpoch(0);

  // Initialize Sensor
  Wire.begin();
  accelSensor.begin(BMA250_range_2g, BMA250_update_time_64ms);

  // Set alarm for initial sleep - allows for delay in recording
  rtc.setAlarmEpoch(rtc.getEpoch() + INITIAL_SLEEP_TIME * SECONDS_IN_HOUR);
  rtc.enableAlarm(rtc.MATCH_YYMMDDHHMMSS);
  LowPower.sleep();

  // Record time at wakeup, record initial timestamp time
  wakeupTime = rtc.getEpoch();
  prevTime = rtc.getEpoch() - EPOCH_OFFSET;
}
```

The setup section configures the SD card, initializes the real-time clock (RTC), and sets up the accelerometer. Lastly, it establishes the initial sleep time, allowing for a delayed start to data collection while keeping the data logger in a very low power sleep mode. 

The data logger enters sleep mode when it reaches the `LowPower.sleep();` line and wakes up after the alarm period. Upon waking, it records the time from the real-time clock to track timestamps for the collected data.

Here’s the complete content formatted as a Markdown file:

#### Loop

```cpp
void loop() {

  // Data array
  double data[DATA_COUNT];

  // Get data
  getAccel(data);

  // Create timestamp
  String timeStampNow = getTimestamp(&prevTime);

  // Record timestamp and data
  writeData(timeStampNow, data);

  // Sleep between samples
  Watchdog.sleep(SLEEP_BETWEEN_SAMPLES);

  int currentTime = rtc.getEpoch();
  // If current recording time interval is over, go to sleep for sleep period
  if (currentTime - wakeupTime > RECORD_INTERVAL) {

    // If data has not been saved to SD card recently, save the data
    if (currentTime - lastSave > SAVE_INTERVAL) {
      lastSave = currentTime;
      file.sync();
    }
    
    rtc.setAlarmEpoch(currentTime + SLEEP_BETWEEN_INTERVAL); // Set alarm SLEEP_BETWEEN_INTERVAL seconds from now
    rtc.enableAlarm(rtc.MATCH_YYMMDDHHMMSS);
    LowPower.sleep();
    wakeupTime = rtc.getEpoch();
  }
}
```

The loop process is as follows: 1) Collect one sample from the accelerometer. 2) Record the data to the microSD card. 3) Sleep before taking the next sample. 4) Upon waking up, check the time. If the recording interval is not over, collect data again. 5) If the recording interval is over, check if it is time to save the data permanently to the SD card ("sync"). This requires a lot of power, so we only do it periodically. Then, set an alarm to sleep until the next recording interval begins. 

## Credits
This project was developed by Brian Katona, with contributions from Wing Sing Chan and Yik-Hei Sung.