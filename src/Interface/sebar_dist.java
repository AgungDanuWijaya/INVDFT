package Interface;

import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.JSchException;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpException;
import java.io.IOException;

/**
 * A program demonstrates how to upload files from local computer to a remote
 * FTP server using Apache Commons Net API.
 *
 * @author www.codejava.net
 */
public class sebar_dist {

    public  void main() throws IOException, JSchException, SftpException {
      
       String Server[] = {"10.42.0.101"};
        String key[] = {"/home/agung/known"};
       
        for (int i = 0; i < Server.length; i++) {
         JSch jsch = new JSch();
        jsch.setKnownHosts(key[i]);
        Session jschSession = jsch.getSession("agung", Server[i]);
        jschSession.setPassword("yut28092018DAM^");
        jschSession.connect();
        ChannelSftp channelSftp = (ChannelSftp) jschSession.openChannel("sftp");
        channelSftp.connect();

        String remoteFile = "/home/agung/incdft/Quantum (copy)/JQC_data/ann_new_dft";
        String localDir = "/home/agung/project/Quantum-20211024T084747Z-001 (2)/Quantum (copy)/JQC_data/ann_new_dft";

        channelSftp.put(localDir, remoteFile );
        System.out.println("main.FTPDownloadFileDemo.main()");
        channelSftp.exit();
        channelSftp.disconnect();
        jschSession.disconnect();	   
        }
        
    }
}