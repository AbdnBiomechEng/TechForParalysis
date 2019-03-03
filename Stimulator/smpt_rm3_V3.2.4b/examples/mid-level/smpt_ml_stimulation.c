#include "smpt_ml_client.h"

static void mid_level_stimulation(const char *port_name);
static void fill_ml_init(Smpt_device * const device, Smpt_ml_init *const ml_init);
static void fill_ml_update(Smpt_device * const device, Smpt_ml_update *const ml_update);
static void fill_ml_get_current_data(Smpt_device * const device, Smpt_ml_get_current_data *const ml_get_current_data);

int main()
{
    /* EDIT: Change to the virtual com port of your device */
    const char *port_name = "COM13";

    mid_level_stimulation(port_name);
    return 0;
}

void mid_level_stimulation(const char *port_name)
{
    Smpt_device device = {0};
    smpt_open_serial_port(&device, port_name);

    Smpt_ml_init ml_init = {0};           /* Struct for ml_init command */
    fill_ml_init(&device, &ml_init);
    smpt_send_ml_init(&device, &ml_init); /* Send the ml_init command to the stimulation unit */

    Smpt_ml_update ml_update = {0};       /* Struct for ml_update command */
    fill_ml_update(&device, &ml_update);
    smpt_send_ml_update(&device, &ml_update);

    Smpt_ml_get_current_data ml_get_current_data = {0};
    fill_ml_get_current_data(&device, &ml_get_current_data);
    smpt_send_ml_get_current_data(&device, &ml_get_current_data);

    smpt_send_ml_stop(&device, smpt_packet_number_generator_next(&device));

    smpt_close_serial_port(&device);
}

void fill_ml_init(Smpt_device *const device, Smpt_ml_init *const ml_init)
{
    /* Clear ml_init struct and set the data */
    smpt_clear_ml_init(ml_init);
    ml_init->packet_number = smpt_packet_number_generator_next(device);
}

void fill_ml_update(Smpt_device *const device, Smpt_ml_update *const ml_update)
{
    /* Clear ml_update and set the data */
    smpt_clear_ml_update(ml_update);
    ml_update->enable_channel[Smpt_Channel_Red] = true;  /* Enable channel red */
    ml_update->packet_number = smpt_packet_number_generator_next(device);

    ml_update->channel_config[Smpt_Channel_Red].number_of_points = 3;  /* Set the number of points */
    ml_update->channel_config[Smpt_Channel_Red].ramp = 3;              /* Three lower pre-pulses   */
    ml_update->channel_config[Smpt_Channel_Red].period = 20;           /* Frequency: 50 Hz */

    /* Set the stimulation pulse */
    /* First point, current: 20 mA, positive, pulse width: 200 µs */
    ml_update->channel_config[Smpt_Channel_Red].points[0].current = 20;
    ml_update->channel_config[Smpt_Channel_Red].points[0].time = 200;

    /* Second point, pause 100 µs */
    ml_update->channel_config[Smpt_Channel_Red].points[1].time = 100;

    /* Third point, current: -20 mA, negative, pulse width: 200 µs */
    ml_update->channel_config[Smpt_Channel_Red].points[2].current = -20;
    ml_update->channel_config[Smpt_Channel_Red].points[2].time = 200;
}

void fill_ml_get_current_data(Smpt_device *const device, Smpt_ml_get_current_data *const ml_get_current_data)
{
    ml_get_current_data->packet_number = smpt_packet_number_generator_next(device);
    ml_get_current_data->data_selection[Smpt_Ml_Data_Stimulation] = true; /* get stimulation data */
}
