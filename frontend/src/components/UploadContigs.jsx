import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Options from './Options.jsx';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import SearchBar from './SearchBar.jsx';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import green from '@material-ui/core/colors/green';
// import { Scrollbars } from 'react-custom-scrollbars';
import download from 'downloadjs';


const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(green[600]),
        backgroundColor: green[500],
        '&:hover': {
            backgroundColor:green[600],
        },
    }
})

class Upload_contigs extends React.Component {

    constructor(props) {

        super(props);

        fetch('api/profiling/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };

        // TODO: poor performnce issue with everytime load the component will
        // fetch API once.

        window.databaseName = "";
        window.fileName = [];

        this.state = {
            to: "/",
            upload_confirm: false,
            switch: false,
        };

        this.djsConfig = {
            dictDefaultMessage:"Drop files or click to upload contigs",
            addRemoveLinks: true,
            maxFilesize:10,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            parallelUploads: 200,
            init:function(){
                this.on("addedfile", function(on_load_header_data){
                    var fileSizeElement = on_load_header_data.previewElement.querySelector(".dz-size");
                    fileSizeElement.parentNode.removeChild(fileSizeElement);
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.batchid);
                });
                this.on("success", function(file){
                    window.fileName.push(file.name);
                    file._removeLink.remove();
                    delete file._removeLink;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.fasta','.fna','.fa'],
            showFiletypeIcon: true,
            postUrl: 'api/profiling/sequence/'
        };

        this.dropzone = null;
    }

    handlePost() {
        
        let fileCheck = this.dropzone.files.length;

        if(fileCheck < 5){
            alert('Please upload at least 5 files');
            return ;
        }

        if(window.databaseName == ""){
            alert('Please choose a database !');
            return ;
        }

        this.dropzone.processQueue();
        this.setState(state => ({ upload_confirm: true , to: '/profile_view' ,
            switch: true }));
    }

    submit(){

        if (this.state.upload_confirm == false){
            alert('Please upload files first! (At least 5 files)');
            return ;
        }

        var scheme = {};
        scheme.occurrence = "95";
        scheme.database = window.databaseName;
        scheme.id = window.batchid;
        fetch('api/profiling/profiling-tree/', {
            method:'POST',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

        window.tabSwitch = true;

    }

    remove(){
        
        this.dropzone.removeAllFiles();
        this.setState(state => ({ to: '/', upload_confirm: false, switch: false }));
        fetch('api/profiling/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };
    }

    query(){

        fetch('api/profiling/profile/' + this.state.querybyID, { method:'GET'})
        .then(response => response.json())
        .catch(error => alert("Data not found"));

        fetch('api/profiling/profile/' + this.state.querybyID, { method:'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                profile_result_all: result.file,
                profile_result_zip: result.zip,})));

        fetch('api/dendrogram/dendrogram/' + this.state.querybyID, { method: 'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                png_file: result.png_file, 
                pdf_file: result.pdf_file,
                svg_file: result.svg_file, 
                emf_file: result.emf_file, 
                newick_file: result.newick_file, })));
    }

    upload_example_data(){

        let mockfile = [
            { name:"S.Agona_01.fa" },
            { name:"S.Agona_02.fa" },
            { name:"S.Agona_03.fa" },
            { name:"S.Enteritidis_01.fa" },
            { name:"S.Enteritidis_02.fa" },
            { name:"S.Enteritidis_03.fa" },
            { name:"S.GoldCoast_01.fa" },
            { name:"S.GoldCoast_02.fa" },
            { name:"S.GoldCoast_03.fa" },
            { name:"S.GoldCoast_04.fa" },
        ];

        let i = 0;
        for(i; i < mockfile.length; i++){
            this.dropzone.emit("addedfile",mockfile[i]);
            this.dropzone.emit("success",mockfile[i]);
            this.dropzone.emit("complete",mockfile[i]);
            this.dropzone.files.push(mockfile[i]);
        };

        let encodeExampleData = [
            require('./static/Example_data/S.Agona_01.fa'), 
            require('./static/Example_data/S.Agona_02.fa'), 
            require('./static/Example_data/S.Agona_03.fa'), 
            require('./static/Example_data/S.Enteritidis_01.fa'), 
            require('./static/Example_data/S.Enteritidis_02.fa'), 
            require('./static/Example_data/S.Enteritidis_03.fa'), 
            require('./static/Example_data/S.Enteritidis_04.fa'), 
            require('./static/Example_data/S.GoldCoast_01.fa'), 
            require('./static/Example_data/S.GoldCoast_02.fa'), 
            require('./static/Example_data/S.GoldCoast_03.fa'), 
        ];

        let decodeExampleData = [];
        let j = 0;

        for(j; j < encodeExampleData.length; j++){
            let tmp = encodeExampleData[j].substring(13,encodeExampleData[j].length);
            tmp = window.atob(tmp);
            decodeExampleData.push(tmp);
        };

        let agona01 = new File([decodeExampleData[0]],'S.Agona_01.fa');
        let agona02 = new File([decodeExampleData[1]],'S.Agona_02.fa');
        let agona03 = new File([decodeExampleData[2]],'S.Agona_03.fa');
        let enteritidis01 = new File([decodeExampleData[3]],'S.Enteritidis_01.fa');
        let enteritidis02 = new File([decodeExampleData[4]],'S.Enteritidis_02.fa');
        let enteritidis03 = new File([decodeExampleData[5]],'S.Enteritidis_03.fa');
        let enteritidis04 = new File([decodeExampleData[6]],'S.Enteritidis_04.fa');
        let goldCoast01 = new File([decodeExampleData[7]],'S.GoldCoast_01.fa');
        let goldCoast02 = new File([decodeExampleData[8]],'S.GoldCoast_02.fa');
        let goldCoast03 = new File([decodeExampleData[9]],'S.GoldCoast_03.fa');

        let decodeExampleFile = [ agona01, agona02, agona03, enteritidis01, enteritidis02,
            enteritidis03, enteritidis04, goldCoast01, goldCoast02, goldCoast03 ];

        let k = 0;
        let upload_status = [];

        for(k; k < decodeExampleFile.length; k++){
            let form = new FormData();
            form.append('file',decodeExampleFile[k]);
            form.append('batch_id',window.batchid);

            fetch('api/profiling/sequence/', {
                method:'POST',
                body:form ,
            }).then(res => upload_status.push(res.status));
        };

        window.databaseName = "Salmonella_enterica";
        window.tabSwitch = true;

        this.setState(state => ({ switch: true }));

        function trigger_celery(){

            if( upload_status.length == 10 ){
                let scheme = {};
                scheme.occurrence = "95";
                scheme.database = window.databaseName;
                scheme.id = window.batchid;
                fetch('api/profiling/profiling-tree/', {
                    method:'POST',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });
                clearInterval(interval);
            }
        };

        let interval = setInterval(trigger_celery,1500);
    }

    back(){
        window.location.reload();
    }

    
    render() {

        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        if(this.state.profile_result_all == undefined){
            return (
            <div>
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <SearchBar
                        onChange = {(value) => this.setState({ querybyID: value })}
                        onRequestSearch={this.query.bind(this)}
                        placeholder = {"Input batch ID here to query data"}
                        style = {{
                            width: '90%',
                            margin: '0 auto',
                        }}
                    />
                </div>
                <br />
                <div style={{ width:'97%', display:'flex', justifyContent:'flex-end', 
                alignItems:'flex-end'}}>
                    <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                            Remove all files
                            &nbsp;&nbsp;
                            <DeleteIcon />
                    </Button>
                </div>
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                        djsConfig={djsConfig} />
                </div>
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <font> SUGGESTION: Please use &nbsp;
                        <a href="http://cab.spbu.ru/software/spades/" target="_blank">SPAdes</a>
                        &nbsp; to assembly before upload. 
                    </font>
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Options switch={this.state.switch} />
                </div>
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" className ={classes.cssRoot} 
                     onClick={this.handlePost.bind(this)}>
                        Upload
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Link to={this.state.to} style={{ textDecoration:'none' }}>
                        <Button variant="contained" color="primary" onClick={this.submit.bind(this)}>
                            profiling
                            &nbsp;&nbsp;
                            <Icon>send</Icon>
                        </Button>
                    </Link>
                </div>
                <br />
                <br />
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" color="default" 
                     onClick={this.upload_example_data.bind(this)}>
                        <Link to="/profile_view" style={{ textDecoration:'none', color:'#000' }}>
                            upload sample data
                            &nbsp;&nbsp;
                            <CloudUploadIcon />
                        </Link>
                    </Button>
                </div>
                <br />
                <br />
                <br />
                <br />
            </div>
            );
        }else{
            return (
                    <div>
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <a download href={this.state.profile_result_all} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.tsv)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.profile_result_zip} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.zip)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <img src={this.state.svg_file} />
                        </div>
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <font size="4">Download Dendrogram</font>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            
                            <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Png 
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Pdf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Svg
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                emf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                newick
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                            <Link to="/" style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default" onClick={this.back.bind(this)}>
                                    <ReplyIcon />
                                    &nbsp;&nbsp;
                                    Back
                                </Button>
                            </Link>
                        </div>
                        <br />
                        <br />
                    </div>
                );
        }
        
    }
}


export default withStyles(styles)(Upload_contigs);