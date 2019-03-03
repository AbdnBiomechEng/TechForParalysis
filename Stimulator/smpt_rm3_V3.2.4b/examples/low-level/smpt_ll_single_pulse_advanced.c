#include "smpt_ll_client.h"

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>


static void init_stimulation(Smpt_device *const device);
static void get_battery_status(Smpt_device *const device);
static void single_pulse(Smpt_device *const device);
static void stop_stimulation(Smpt_device *const device);

static void start_timing(void);
static int  stop_timing(void);

struct timeval start_, stop_;
long u1_ = 0, u2_ = 0;

int main(void)
{
	/* EDIT: Change to the virtual com port of your device */
    const char *port_name = "COM5";

    Smpt_device device = {0};
    smpt_open_serial_port(&device, port_name);

    init_stimulation(  &device);
    get_battery_status(&device);
    single_pulse(      &device);    /* Sends one bi-phasic impulse and checks the response */
    stop_stimulation(  &device);

    smpt_close_serial_port(&device);

	return 0;
}

void init_stimulation(Smpt_device *const device)
{
    Smpt_ll_init     ll_init     = {0};  /* Struct for ll_init command */
	Smpt_ll_init_ack ll_init_ack = {0};  /* Struct for ll_init_ack response */
    Smpt_ack ack                 = {0};  /* Struct for general response */

	smpt_clear_ll_init(&ll_init);
    ll_init.packet_number = smpt_packet_number_generator_next(device);

    printf("SMPT init_stimulaton(): send Ll_init command ...\n");
    start_timing();
    smpt_send_ll_init(device, &ll_init);   /* Send the ll_init command to the stimulation unit */

    while (!smpt_new_packet_received(device)) { /* busy waits for Ll_init_ack response */}
    printf("SMPT init_stimulaton(): Ll_init_ack received, took %d ms...\n", stop_timing());

	smpt_clear_ack(&ack);
    smpt_last_ack(device, &ack);
	if (ack.command_number == Smpt_Cmd_Ll_Init_Ack)
	{
        smpt_get_ll_init_ack(device, &ll_init_ack);  /* Writes the received data into ll_init_ack */
		if ( ll_init_ack.result == Smpt_Result_Successful )
		{
            printf("SMPT init_stimulation(): Ll_init command was successful\n");
		}
		else
		{
            printf("SMPT init_stimulation(): Ll_init command failed! Return code: %d\n", (int)ll_init_ack.result);
		}
	}
	else
	{
        printf("SMPT init_stimulation(): Unexpected ack received! Command ID: %d\n", (int)ack.command_number);
	}
	printf("\n\n");
}

void get_battery_status(Smpt_device *const device)
{
    Smpt_ll_get_status_ack ll_get_status_ack = {0}; /* Struct for ll_get_status_ack response */
    Smpt_ack ack = {0};                             /* Struct for general response */

    printf("SMPT get_battery_status(): send Get_Status command ...\n");
    start_timing();
    smpt_send_ll_get_status(device, smpt_packet_number_generator_next(device));

    while (!smpt_new_packet_received(device)) { /* busy waits for Ll_status_ack response */}
    printf("SMPT get_battery_status(): Status ack received, took %d ms...\n", stop_timing());

    smpt_clear_ack(&ack);
    smpt_last_ack(device, &ack);
    if (ack.command_number == Smpt_Cmd_Ll_Get_Status_Ack)
    {
        smpt_get_ll_get_status_ack(device, &ll_get_status_ack);  /* Writes the received data into ll_init_ack */

        if ( ll_get_status_ack.result == Smpt_Result_Successful )
        {
            printf("SMPT get_battery_status(): Get_Status command was successful:\n");
            printf("    Battery Voltage: %fV\n    Battery Level: %d%%\n",
                   (double) (ll_get_status_ack.battery_voltage /1000.0),
                   (int) ll_get_status_ack.battery_level);
        }
        else
        {
            printf("SMPT get_battery_status(): Get_Status command failed! Return code: %d\n", (int)ll_get_status_ack.result);
        }
    }
    else
    {
        printf("SMPT get_battery_status(): Unexpected ack received! Command ID: %d\n", (int)ack.command_number);
    }
    printf("\n\n");
}


void single_pulse(Smpt_device *const device)
{
    Smpt_ll_channel_config ll_channel_config         = {0}; /* Struct for ll_channel_config command */
    Smpt_ll_channel_config_ack ll_channel_config_ack = {0}; /* Struct for the ll_channel_config_ack response */
    Smpt_ack ack = {0};                                     /* Struct for general response */

    printf("SMPT single_pulse(): send Puls_Config command ...\n");
    start_timing();

	/* Set the data */
	ll_channel_config.enable_stimulation = true;
	ll_channel_config.channel = Smpt_Channel_Red;  /* Use red channel */
	ll_channel_config.number_of_points = 3;         /* Set the number of points*/
    ll_channel_config.packet_number = smpt_packet_number_generator_next(device);

	/* Set the stimulation pulse */
	/* First point, current: 20 mA, positive, pulse width: 200 µs */
	ll_channel_config.points[0].current =  20;
	ll_channel_config.points[0].time    = 200;

	/* Second point, pause 100 µs */
	ll_channel_config.points[1].time = 100;

	/* Third point, current: -20 mA, negative, pulse width: 200 µs */
	ll_channel_config.points[2].current = -20;
	ll_channel_config.points[2].time    = 200;

	/* Send the ll_channel_list command to the stimulation unit */
    smpt_send_ll_channel_config(device, &ll_channel_config);

    while (!smpt_new_packet_received(device)) { /* busy waits for ll_stop_ack */ }
	/* Get start time for timing measurements */
    start_timing();

	smpt_clear_ack(&ack);
    smpt_last_ack(device, &ack);
	if (ack.command_number == Smpt_Cmd_Ll_Channel_Config_Ack)
	{
        smpt_get_ll_channel_config_ack(device, &ll_channel_config_ack);  /* Writes the received data into ll_channel_config_ack */
		/* Give Feedback about the initialisation */
		if ( ll_channel_config_ack.result == Smpt_Result_Successful )
		{
            printf("SMPT single_pulse_test: Puls_Config command was successful! (took %fms)\n", stop_timing());
		}
		else
		{
			printf("SMPT single_pulse_test: Puls_Config command failed! Return code: %d; Electrode Error Code: %d\n", (int)ll_channel_config_ack.result, (int)ll_channel_config_ack.electrode_error);
            stop_timing();
		}
	}
	else
	{
		printf("SMPT single_pulse_test: Unexpected ack received! Command ID: %d\n", (int)ack.command_number);
	}
	printf("\n\n");
}

void stop_stimulation(Smpt_device *const device)
{
    Smpt_ack ack = {0};

    printf("SMPT wait_for_ACK_test: send Stop command ...\n");
    start_timing();
    smpt_send_ll_stop(device, smpt_packet_number_generator_next(device));

    while (!smpt_new_packet_received(device)) { /* busy waits for Ll_stop_ack response */}

    smpt_clear_ack(&ack);
    smpt_last_ack(device, &ack);
    if (ack.command_number == Smpt_Cmd_Ll_Stop_Ack)
    {
        /* Ll_stop_ack has been received */
        printf("SMPT wait_for_ACK_test: Stop ack received, took %d ms...\n", stop_timing());
    }
    else
    {
        printf("SMPT wait_for_ACK_test: Unexpected ack received! Command ID: %d\n", (int)ack.command_number);
        stop_timing();
    }

    printf("\n\n");
}

void start_timing(void)
{
    /* Get start time for timing measurements */
    gettimeofday(&start_, NULL);
    u1_ = start_.tv_sec * 1000 + start_.tv_usec / 1000;
}

int stop_timing(void)
{
    /* Get start time for timing measurements */
    gettimeofday(&stop_, NULL);
    u2_ = stop_.tv_sec*1000 + stop_.tv_usec/1000;

    return (int)(u2_-u1_);
}
