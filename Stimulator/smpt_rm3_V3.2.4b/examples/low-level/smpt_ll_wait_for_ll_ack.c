#include "smpt_ll_client.h"

static void smpt_wait_for_ll_ack(const char *port_name);

int main(void)
{
    /* EDIT: Change to the virtual com port of your device */
    const char *port_name = "COM5";

    smpt_wait_for_ll_ack(port_name);
    return 0;
}

void smpt_wait_for_ll_ack(const char *port_name)
{
    uint8_t packet_number = 0;  /* The packet_number can be used for debugging purposes */
    Smpt_ll_init ll_init = {0};          /* Struct for ll_init command */
    Smpt_ll_init_ack ll_init_ack = {0};  /* Struct for ll_init_ack response */
    Smpt_ack ack = {0};            /* Struct for general response */

    Smpt_device device = {0};
    smpt_open_serial_port(&device, port_name);

    /* Clear ll_init struct and set the data */
    smpt_clear_ll_init(&ll_init);
    ll_init.packet_number = packet_number;

    /* Send the ll_init command to the stimulation unit */
    smpt_send_ll_init(&device, &ll_init);

    packet_number++;

    while (!smpt_new_packet_received(&device)) { /* busy waits for Ll_init_ack response */}

    smpt_clear_ack(&ack);
    smpt_last_ack(&device, &ack);
    if (ack.command_number == Smpt_Cmd_Ll_Init_Ack)
    {
        smpt_get_ll_init_ack(&device, &ll_init_ack);  /* Writes the received data into ll_init_ack */
    }

    /* Send the Ll_stop command to stimulation unit */
    smpt_send_ll_stop(&device, packet_number);

    while (!smpt_new_packet_received(&device)) { /* busy waits for ll_stop_ack */ }

    smpt_clear_ack(&ack);
    smpt_last_ack(&device, &ack);
    if (ack.command_number == Smpt_Cmd_Ll_Stop_Ack)
    {
        /* Ll_stop_ack has been received */
    }

    smpt_close_serial_port(&device);
}
